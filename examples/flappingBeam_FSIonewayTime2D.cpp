/// This is a oneway FSI coupling benchmark based on the fluid-structure interaction benchmark FSI2 from this paper:
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
              index_t flowIter, index_t beamIter)
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

    // print: simTime drag lift pressureDiff dispAx dispAy aleNorm aleTime flowTime beamTime flowIter beamIter
    ofs << simTime << " " << force.at(0) << " " << force.at(1) <<  " " << presFront.at(0)-presBack.at(0) << " "
        << dispA.at(0) << " " << dispA.at(1) << " " << normL2(geoALE,dispALE) << " "
        << aleTime << " " << flowTime << " " << beamTime << " " << flowIter << " " << beamIter << std::endl;
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the one-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // beam parameters
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t densitySolid = 1.0e3;
    real_t loading = 2;
    // flow parameters
    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flow.xml";
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    // ALE parameters
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.5; // local stiffening for ALE
    index_t ALEmethod = ale_method::TINE;
    bool check = true;
    // space discretization
    index_t numUniRef = 3;
    // time integration
    real_t timeStep = 0.01;
    real_t timeSpan = 3.;
    real_t thetaFluid = 0.5;
    bool imexOrNewton = false;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 1000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the one-way fluid-structure interaction solver in 2D.");
    cmd.addReal("l","load","Gravitation acceleration acting on the beam",loading);
    cmd.addReal("m","meanvelocity","Average inflow velocity",meanVelocity);
    cmd.addReal("v","viscosity","Fluid kinematic viscosity",viscosity);
    cmd.addReal("x","chi","Local stiffening degree for ALE",meshStiff);
    cmd.addSwitch("c","check","Check bijectivity of the ALE displacement field",check);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("f","thetafluid","One-step theta scheme for the fluid: 0 - exp.Euler, 1 - imp.Euler, 0.5 - Crank-Nicolson",thetaFluid);
    cmd.addSwitch("i","intergration","Time integration scheme for the fluid: false = IMEX (default), true = Newton",imexOrNewton);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of sampling points per patch for Paraview (0 = no plotting)",numPlotPoints);
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

    // creating matching bases
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
    gsConstantFunction<> gFlow(0.,0.,2);
    gsConstantFunction<> gBeam(0.,loading*densitySolid,2);
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
    // FSI sides: flow first, then beam
    gsInterfaceFSI interface;
    interface.addSide(0,boundary::south,0,boundary::north);
    interface.addSide(1,boundary::north,0,boundary::south);
    interface.addSide(2,boundary::west,0,boundary::east);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,gFlow);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,gFlow);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler,&velALE,&patchesALE);
    nsTimeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gFlow);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
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

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velFlow);
    gsField<> pressureField(nsAssembler.patches(),presFlow);
    gsField<> displacementField(geoBeam,dispBeam);
    gsField<> aleField(geoALE,dispALE);

    // creating containers to plot fields to Paraview files
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSIow_flow");
    gsParaviewCollection collectionBeam("flappingBeam_FSIow_beam");
    gsParaviewCollection collectionALE("flappingBeam_FSIow_ALE");

    std::ofstream logFile;
    logFile.open("flappingBeam_FSIow.txt");
    logFile << "# simTime drag lift pressureDiff dispAx dispAy aleNorm aleTime flowTime beamTime flowIter beamIter\n";

    gsProgressBar bar;
    gsStopwatch iterClock, totalClock;

    //=============================================//
                   // Initial conditions //
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
    moduleALE.constructSolution(dispALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSIow_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSIow_beam",collectionBeam,0,numPlotPoints);
        //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSIow_ALE",collectionALE,0,numPlotPoints);
        plotDeformation(geoALE,dispALE,"flappingBeam_FSIow_ALE",collectionALE,0);
    }
    writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,0.,0.,0.,0.,0,0);

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

        // BEAM
        iterClock.restart();
        gsMultiPatch<> dispDiff;
        elAssembler.constructSolution(elTimeSolver.displacementVector(),dispDiff);
        elTimeSolver.makeTimeStep(tStep);
        elAssembler.constructSolution(elTimeSolver.displacementVector(),dispBeam);
        dispDiff.patch(0).coefs() = dispBeam.patch(0).coefs() - dispDiff.patch(0).coefs();
        timeBeam += iterClock.stop();

        // ALE
        iterClock.restart();
        index_t badPatch = moduleALE.updateMesh();
        moduleALE.constructSolution(velALE);
        for (index_t p = 0; p < velALE.nPatches(); ++p)
        {
            // construct ALE difference
            velALE.patch(p).coefs() -= dispALE.patch(p).coefs();
            // update flow geo
            nsAssembler.patches().patch(p+3).coefs() += velALE.patch(p).coefs();
            nsMassAssembler.patches().patch(p+3).coefs() += velALE.patch(p).coefs();
            // construct ALE velocity
            velALE.patch(p).coefs() /= tStep;
        }
        // construct ALE displacement
        moduleALE.constructSolution(dispALE);
        timeALE += iterClock.stop();
        // test if any patch is not bijective
        if (badPatch != -1)
        {
            gsInfo << "\n Bad patch: " << badPatch << std::endl;
            break;
        }

        // FLOW
        iterClock.restart();
        nsAssembler.setFixedDofs(3,boundary::south,velALE.patch(0).boundary(boundary::south)->coefs());
        nsAssembler.setFixedDofs(4,boundary::north,velALE.patch(1).boundary(boundary::north)->coefs());
        nsAssembler.setFixedDofs(5,boundary::west,velALE.patch(2).boundary(boundary::west)->coefs());
        nsTimeSolver.makeTimeStep(tStep);
        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velFlow,presFlow);
        timeFlow += iterClock.stop();

        // Iteration end
        simTime += tStep;
        numTimeStep++;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSIow_flow",collectionFlow,numTimeStep,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSIow_beam",collectionBeam,numTimeStep,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSIow_ALE",collectionALE,numTimeStep,numPlotPoints);
            plotDeformation(geoALE,dispALE,"flappingBeam_FSIow_ALE",collectionALE,numTimeStep);
        }
        writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,
                 simTime,timeALE,timeFlow,timeBeam,
                 nsTimeSolver.numberIterations(),elTimeSolver.numberIterations());
    }

    //=============================================//
                   // Final touches //
    //=============================================//

    gsInfo << "Complete in: " << secToHMS(totalClock.stop())
           << ", ALE time: " << secToHMS(timeALE)
           << ", flow time: " << secToHMS(timeFlow)
           << ", beam time: " << secToHMS(timeBeam) << std::endl;

    if (numPlotPoints > 0)
    {
        collectionFlow.save();
        collectionBeam.save();
        collectionALE.save();
        gsInfo << "Open \"flappingBeam_FSIow_*.pvd\" in Paraview for visualization.\n";
    }
    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_FSIow.txt\".\n";

    return 0;
}
