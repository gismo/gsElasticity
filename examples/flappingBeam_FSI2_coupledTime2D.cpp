/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsALE.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

void writeLog(std::ofstream & ofs, const gsNsAssembler<real_t> & assemblerFlow,
              const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure,
              const gsMultiPatch<> & displacementBeam, const gsMultiPatch<> & geoALE, const gsMultiPatch<> & dispALE,
              real_t simTime, real_t aleTime, real_t flowTime, real_t beamTime,
              index_t couplingIter, index_t flowIter, index_t beamIter,
              real_t omega, real_t resAbs, real_t resRel)
{
    // computing force acting on the surface of the submerged structure
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assemblerFlow.computeForce(velocity,pressure,bdrySides);

    // compute the pressure difference between the front and the end points of the structure
    gsMatrix<> presFront(2,1);
    presFront << 1.,0.5;
    presFront = pressure.patch(0).eval(presFront);
    gsMatrix<> presBack(2,1);
    presBack << 0.,0.5;
    presBack = pressure.patch(5).eval(presBack);

    // compute displacement of the beam point A
    gsMatrix<> dispA(2,1);
    dispA << 1.,0.5;
    dispA = displacementBeam.patch(0).eval(dispA);

    // print: 1-simTime 2-drag 3-lift 4-pressureDiff 5-dispAx 6-dispAy 7-aleNorm
    //        8-aleTime 9-flowTime 10-beamTime 11-couplingIter 12-flowIter 13-beamIter
    //        14-omega 15-resAbs 16-resRel
    ofs << simTime << " " << force.at(0) << " " << force.at(1) << " " << presFront.at(0)-presBack.at(0) << " "
        << dispA.at(0) << " " << dispA.at(1) << " " << normL2(geoALE,dispALE) << " "
        << aleTime << " " << flowTime << " " << beamTime << " "
        << couplingIter << " " << flowIter << " " << beamIter << " "
        << omega << " " << resAbs << " " << resRel << std::endl;
}

void formVector(const gsMultiPatch<> & disp, gsMatrix<> & vector)
{
    gsMatrix<> north = disp.patch(0).boundary(boundary::north)->coefs();
    gsMatrix<> south = disp.patch(0).boundary(boundary::south)->coefs();
    gsMatrix<> east = disp.patch(0).boundary(boundary::east)->coefs();
    index_t northSize = north.rows();
    index_t southSize = south.rows();
    index_t eastSize = east.rows();
    index_t totalSize = northSize+southSize+eastSize;
    index_t dim = north.cols();

    vector.setZero(totalSize*dim,1);
    for (index_t d = 0; d < dim;++d)
    {
        vector.middleRows(totalSize*d,northSize) = north.col(d);
        vector.middleRows(totalSize*d+northSize,southSize) = south.col(d);
        vector.middleRows(totalSize*d+northSize+southSize,eastSize) = east.col(d);
    }
}


real_t aitken(gsMultiPatch<> & dispA, gsMultiPatch<> & dispB,gsMultiPatch<> & dispB2, gsMultiPatch<> & dispC, real_t & omega,bool & converged, real_t & resNorInit)
{
    gsMatrix<> vecA,vecB,vecB2,vecC;
    formVector(dispA,vecA);
    formVector(dispB,vecB);
    formVector(dispB2,vecB2);
    formVector(dispC,vecC);


    omega = -1*omega*((vecB2-vecA).transpose()*(vecC-vecB -vecB2+vecA))(0,0)/((vecC-vecB -vecB2+vecA).transpose()*(vecC-vecB -vecB2+vecA))(0,0);
    dispA.patch(0).coefs() = dispB.patch(0).coefs();
    dispB.patch(0).coefs() += omega*(dispC.patch(0).coefs()-dispB.patch(0).coefs());
    dispB2.patch(0).coefs() = dispC.patch(0).coefs();

    real_t resnor = ((vecC-vecB)*omega).norm();
    //gsInfo <<"Resnor abs " << resnor << ", resnor rel " << resnor/resNorInit << std::endl;
    if (resnor < 1e-10 || resnor/resNorInit < 1e-6)
        converged = true;

    //gsInfo << "Omega " << omega << std::endl;
    return resnor;
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the two-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // beam parameters
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t densitySolid = 1.0e4;
    // flow parameters
    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flow.xml";
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    // ALE parameters
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.5;
    index_t ALEmethod = ale_method::TINE;
    bool check = true;
    // space discretization
    index_t numUniRef = 3;
    // time integration
    real_t timeStep = 0.01;
    real_t timeSpan = 15.;
    real_t thetaFluid = 0.5;
    real_t thetaSolid = 1.;
    index_t maxCouplingIter = 10;
    bool imexOrNewton = false;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 0.;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("x","chi","Local stiffening degree for ALE",meshStiff);
    cmd.addSwitch("c","check","Check bijectivity of the ALE displacement field",check);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addInt("i","iter","Number of coupling iterations",maxCouplingIter);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);
    // I deform only three flow patches adjacent to the FSI interface
    gsMultiPatch<> geoALE;
    for (index_t p = 0; p < 3; ++p)
        geoALE.addPatch(geoFlow.patch(p+3).clone());
    geoALE.computeTopology();
    // correspondence mapping between the flow and the ALE patches
    std::vector<std::pair<index_t,index_t> > patchesALE;
    for (index_t p = 0; p < 3; ++p)
        patchesALE.push_back(std::pair<index_t,index_t>(p+3,p));

    // creating bases
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        geoFlow.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisPressure(geoFlow);
    // I use subgrid elements (because degree elevation is not implemented for geometries, so I cant use Taylor-Hood)
    basisDisplacement.uniformRefine();
    geoALE.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velFlow, presFlow, dispBeam, dispALE, velALE;
    // boundary conditions: flow
    gsBoundaryConditions<> bcInfoFlow;
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0);
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfoFlow.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(5,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::north,condition_type::dirichlet,0,d);
    }
    // boundary conditions: beam
    gsBoundaryConditions<> bcInfoBeam;
    for (index_t d = 0; d < 2; ++d)
        bcInfoBeam.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
    gsFsiLoad<real_t> fSouth(geoALE,dispALE,1,boundary::north,
                             velFlow,presFlow,4,viscosity,densityFluid);
    gsFsiLoad<real_t> fEast(geoALE,dispALE,2,boundary::west,
                            velFlow,presFlow,5,viscosity,densityFluid);
    gsFsiLoad<real_t> fNorth(geoALE,dispALE,0,boundary::south,
                             velFlow,presFlow,3,viscosity,densityFluid);
    bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&fSouth);
    bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&fEast);
    bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&fNorth);

    // beam to ALE interface
    gsInterfaceFSI interface;
    interface.addSide(0,boundary::south,0,boundary::north);
    interface.addSide(1,boundary::north,0,boundary::south);
    interface.addSide(2,boundary::west,0,boundary::east);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,g);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,g);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler,&velALE,&patchesALE);
    nsTimeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    nsTimeSolver.options().setReal("Theta",thetaFluid);
    nsTimeSolver.options().setSwitch("ALE",true);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    elTimeSolver.options().setReal("Beta",thetaSolid/2);
    elTimeSolver.options().setReal("Gamma",thetaSolid);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // mesh deformation module
    gsALE<real_t> moduleALE(geoALE,dispBeam,interface,ale_method::method(ALEmethod));
    moduleALE.options().setReal("LocalStiff",meshStiff);
    moduleALE.options().setReal("PoissonsRatio",meshPR);
    moduleALE.options().setSwitch("Check",check);
    gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // beam stress field
    gsPiecewiseFunction<> stresses;

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velFlow);
    gsField<> pressureField(nsAssembler.patches(),presFlow);
    gsField<> displacementField(geoBeam,dispBeam);
    gsField<> stressField(geoBeam,stresses,true);
    gsField<> aleField(geoALE,dispALE);

    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    fieldsBeam["von Mises"] = &stressField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSI2_flow");
    gsParaviewCollection collectionBeam("flappingBeam_FSI2_beam");
    gsParaviewCollection collectionALE("flappingBeam_FSI2_ALE");

    std::ofstream logFile;
    logFile.open("flappingBeam_FSI2.txt");
    logFile << "# simTime drag lift presDiff dispAx dispAy aleNorm aleTime flowTime"
            << " beamTime couplingIter flowIter beamIter omega resAbs resRel\n";

    gsProgressBar bar;
    gsStopwatch totalClock, iterClock;

    //=============================================//
                   // Initial condtions //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.getFixedDofs(0,boundary::west,inflowDDoFs);
    nsAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());

    elTimeSolver.setDisplacementVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));
    elTimeSolver.setVelocityVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));

    // plotting initial condition
    nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velFlow,presFlow);
    elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispBeam);
    elAssembler.constructCauchyStresses(dispBeam,stresses,stress_components::von_mises);
    moduleALE.constructSolution(dispALE);

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,0,numPlotPoints);
        //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,0,numPlotPoints);
        plotDeformation(geoALE,dispALE,"flappingBeam_FSI2_ALE",collectionALE,0);
    }
    writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,0.,0.,0.,0.,0,0,0,1.,0.,0.);

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    real_t simTime = 0.;
    real_t numTimeStep = 0;
    real_t timeALE = 0.;
    real_t timeFlow = 0.;
    real_t timeBeam = 0.;

    totalClock.restart();

    gsInfo << "Running the simulation...\n";
    while (simTime < timeSpan)
    {
        bar.display(simTime/timeSpan);
        // change time step for the initial warm-up phase
        real_t tStep = (warmUp && simTime < 2.) ? 0.1 : timeStep;
        // smoothly change the inflow boundary condition
        if (simTime < 2.)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*(simTime+tStep)/2.))/2);

        index_t couplingIter = 0;
        bool converged = false;

        elTimeSolver.saveState();
        moduleALE.saveState();
        nsTimeSolver.saveState();

        gsMultiPatch<> dispOldOld, dispNew,dispDiff, dispGuess;
        real_t omega = 1.;
        real_t resNorInit = 1.;
        real_t resAbs = 0.;

        while (couplingIter < maxCouplingIter && !converged)
        {
            // BEAM
            iterClock.restart();
            if (couplingIter > 0)
                elTimeSolver.recoverState();
            elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispDiff);
            elTimeSolver.makeTimeStep(tStep);
            if (couplingIter == 0)
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispOldOld);
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispBeam);
            }
            else if (couplingIter == 1)
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispBeam);
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispGuess);
                gsMatrix<> vecA, vecB;
                formVector(dispOldOld,vecA);
                formVector(dispBeam,vecB);
                resNorInit = (vecB-vecA).norm();
            }
            else
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispNew);
                resAbs = aitken(dispOldOld,dispBeam,dispGuess,dispNew,omega,converged,resNorInit);
            }
            //elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
            dispDiff.patch(0).coefs() = dispBeam.patch(0).coefs() - dispDiff.patch(0).coefs();
            timeBeam += iterClock.stop();

            // ALE
            iterClock.restart();
            // undo last ALE
            for (index_t p = 0; p < 3; ++p)
            {
                nsAssembler.patches().patch(p+3).coefs() -= dispALE.patch(p).coefs();
                nsMassAssembler.patches().patch(p+3).coefs() -= dispALE.patch(p).coefs();
            }
            // recover ALE at the start of timestep
            if (couplingIter > 0)
                moduleALE.recoverState();
            moduleALE.constructSolution(dispALE);
            // update
            index_t badPatch = moduleALE.updateMesh();
            // construct ALE velocity
            moduleALE.constructSolution(velALE);
            for (index_t p = 0; p < velALE.nPatches(); ++p)
            {
                velALE.patch(p).coefs() -= dispALE.patch(p).coefs();
                velALE.patch(p).coefs() /= tStep;
            }
            // construct ALE displacement
            moduleALE.constructSolution(dispALE);
            // apply new ALE to the flow domain
            for (index_t p = 0; p < 3; ++p)
            {
                nsAssembler.patches().patch(p+3).coefs() += dispALE.patch(p).coefs();
                nsMassAssembler.patches().patch(p+3).coefs() += dispALE.patch(p).coefs();
            }
            timeALE += iterClock.stop();

            // test if any patch is not bijective
            if (badPatch != -1)
            {
                gsInfo << "\n Bad patch: " << badPatch << std::endl;
                goto labelend;
            }

            // FLOW
            iterClock.restart();
            nsAssembler.setFixedDofs(3,boundary::south,velALE.patch(0).boundary(boundary::south)->coefs());
            nsAssembler.setFixedDofs(4,boundary::north,velALE.patch(1).boundary(boundary::north)->coefs());
            nsAssembler.setFixedDofs(5,boundary::west,velALE.patch(2).boundary(boundary::west)->coefs());
            if (couplingIter > 0)
                nsTimeSolver.recoverState();
            nsTimeSolver.makeTimeStep(tStep);
            nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velFlow,presFlow);
            timeFlow += iterClock.stop();

            ++couplingIter;
        }

        // Iteration end
        simTime += tStep;
        numTimeStep++;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,numTimeStep,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,numTimeStep,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,numTimeStep,numPlotPoints);
            plotDeformation(geoALE,dispALE,"flappingBeam_FSI2_ALE",collectionALE,numTimeStep);
        }
        writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,
                 simTime,timeALE,timeFlow,timeBeam, couplingIter,
                 nsTimeSolver.numberIterations(),elTimeSolver.numberIterations(),
                 omega,resAbs,resAbs/resNorInit);
    }

    //=============================================//
                   // Final touches //
    //=============================================//

    labelend:
    gsInfo << "Complete in: " << secToHMS(totalClock.stop())
           << ", ALE time: " << secToHMS(timeALE)
           << ", flow time: " << secToHMS(timeFlow)
           << ", beam time: " << secToHMS(timeBeam) << std::endl;

    if (numPlotPoints > 0)
    {
        collectionFlow.save();
        collectionBeam.save();
        collectionALE.save();
        gsInfo << "Open \"flappingBeam_FSI2_*.pvd\" in Paraview for visualization.\n";
    }
    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_FSI2.txt\".\n";

    return 0;
}
