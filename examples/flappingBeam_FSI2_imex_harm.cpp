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
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElPoissonAssembler.h>

using namespace gismo;

void validation(std::ofstream & ofs, const gsNsAssembler<real_t> & assembler, real_t time,
                const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure,const gsMultiPatch<> & displacement)
{
    // computing force acting on the surface of the cylinder
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);
    // print time-drag-lift-pressure-YdispA
    gsMatrix<> disp(2,1);
    disp << 1.,0.5;
    disp = displacement.patch(0).eval(disp);
    gsMatrix<> vel(2,1);
    vel << 0.,0.5;
    vel = velocity.patch(5).eval(vel);

    ofs << time << " " << force(0,0) << " " << force(1,0) << " "
        << disp(0,0) << " " << disp(1,0) << " " << vel(0,0) << " " << vel(1,0) << std::endl;
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


void aitken(gsMultiPatch<> & dispA, gsMultiPatch<> & dispB,gsMultiPatch<> & dispB2, gsMultiPatch<> & dispC, real_t & omega,bool & converged, real_t & resNorInit)
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
    gsInfo <<"Resnor abs " << resnor << ", resnor rel " << resnor/resNorInit << std::endl;
    if (resnor < 1e-10 || resnor/resNorInit < 1e-6)
        converged = true;

    gsInfo << "Omega " << omega << std::endl;
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the one-way time-dependent fluid-structure interaction solver in 2D.\n";

    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameFlowPart = ELAST_DATA_DIR"/flappingBeam_flowPart.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    real_t densitySolid = 1.0e4;
    real_t timeStep = 0.01;
    real_t timeSpan = 1.;
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.; // local stiffening for ALE
    index_t numIter = 4;
    index_t type = 0;


    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("j","stiff","Local stiffening for ALE",meshStiff);
    cmd.addInt("i","iter","Num iters",numIter);
    cmd.addInt("y","type","Type",type);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoALE; // this is a part of the flow geometry; we deform only this part to save memory and time
    gsReadFile<>(filenameFlowPart, geoALE);
    geoALE.computeTopology(); // just in case
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

    // creating bases
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        geoFlow.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisPressure(geoFlow);
    // we use subgrid elements (degree elevation does not work for domains)
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
    gsConstantFunction<> gALE(0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velocity, pressure, displacement, ALE;
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
    gsFsiLoad<real_t> fSouth(geoALE,ALE,1,boundary::north,
                             velocity,pressure,4,viscosity,densityFluid);
    gsFsiLoad<real_t> fEast(geoALE,ALE,2,boundary::west,
                            velocity,pressure,5,viscosity,densityFluid);
    gsFsiLoad<real_t> fNorth(geoALE,ALE,0,boundary::south,
                             velocity,pressure,3,viscosity,densityFluid);
    bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&fSouth);
    bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&fEast);
    bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&fNorth);

    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoALE;
    for (auto it = geoALE.bBegin(); it != geoALE.bEnd(); ++it)
            bcInfoALE.addCondition(it->patch,it->side(),condition_type::dirichlet,0,0);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,g);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,g);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler);
    nsTimeSolver.options().setInt("Scheme",time_integration::implicit_linear);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: flow mesh
    gsElasticityAssembler<real_t> baseAssembler(geoALE,basisALE,bcInfoALE,g);
    gsInfo << "Initialized elasticity system for ALE with " << baseAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velocity);
    gsField<> pressureField(nsAssembler.patches(),pressure);
    gsField<> displacementField(geoBeam,displacement);
    gsField<> aleField(geoALE,ALE);

    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSI2_flow");
    gsParaviewCollection collectionBeam("flappingBeam_FSI2_beam");
    gsParaviewCollection collectionALE("flappingBeam_FSI2_ALE");

    std::ofstream file;
    file.open("flappingBeam_FSI2.txt");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
                   // Initial condtions //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.getFixedDofs(0,boundary::west,inflowDDoFs);

    // set all Dirichlet DoFs to zero
    nsAssembler.homogenizeFixedDofs(-1);
    nsMassAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());
    nsTimeSolver.initialize();


    // plotting initial condition
    nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
    elAssembler.constructSolution(gsMatrix<>::Zero(elAssembler.numDofs(),1),displacement);
    baseAssembler.constructSolution(gsMatrix<>::Zero(baseAssembler.numDofs(),1),baseAssembler.allFixedDofs(),ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,0,numPlotPoints);
        //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,0,numPlotPoints);
        plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,0);
    }
    validation(file,nsAssembler,0.,velocity,pressure,displacement);

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    std::vector<std::pair<index_t,index_t> > alePatches;
    alePatches.push_back(std::pair<index_t,index_t>(3,0));
    alePatches.push_back(std::pair<index_t,index_t>(4,1));
    alePatches.push_back(std::pair<index_t,index_t>(5,2));

    elTimeSolver.initialize();

    real_t timeALE = 0.;
    real_t timeFlow = 0.;
    real_t timeBeam = 0.;

    gsInfo << "FSI...\n";
    index_t bad;
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        //bar.display(i+1,index_t(timeSpan/timeStep));
        gsInfo << "===========================Step " << i+1 <<"/"<<index_t(timeSpan/timeStep) << std::endl;

        if (type == 0)
            bad = checkGeometry(nsAssembler.patches());
        else
            bad = baseAssembler.checkSolution(ALE);
        if (bad != -1)
        {
            gsInfo << "\n Bad patch: " << bad << std::endl;
            break;
        }

        if (i*timeStep < 2.)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*timeStep*(i+1)/2.))/2);

        index_t iter = 0;
        bool converged = false;

        elTimeSolver.saveState();
        nsTimeSolver.saveState();

        gsMultiPatch<> dispOldOld, dispNew,dispDiff, dispGuess;
        real_t omega = 0.1;
        real_t resNorInit;
        gsMultiPatch<> xiUpdate, etaUpdate;


        while (iter < numIter && !converged)
        {
            // BEAM
            clock.restart();

            if (iter > 0)
                elTimeSolver.recoverState();
            elAssembler.constructSolution(elTimeSolver.displacementVector(),dispDiff);
            elTimeSolver.makeTimeStep(timeStep);
            if (iter == 0)
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispOldOld);
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
            }
            else if (iter == 1)
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispGuess);
                gsMatrix<> vecA, vecB;
                formVector(dispOldOld,vecA);
                formVector(displacement,vecB);
                resNorInit = (vecB-vecA).norm();
            }
            else
            {
                elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),dispNew);
                aitken(dispOldOld,displacement,dispGuess,dispNew,omega,converged,resNorInit);
            }
            //elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
            dispDiff.patch(0).coefs() = displacement.patch(0).coefs() - dispDiff.patch(0).coefs();
            timeBeam += clock.stop();

            // ALE
            clock.restart();

            // undo last ALE
            if (iter > 0)
                for (index_t p = 0; p < 3; ++p)
                {
                    nsAssembler.patches().patch(p+3).coefs().col(0) -= xiUpdate.patch(p).coefs();
                    nsAssembler.patches().patch(p+3).coefs().col(1) -= etaUpdate.patch(p).coefs();
                    ALE.patch(p).coefs().col(0) -= xiUpdate.patch(p).coefs();
                    ALE.patch(p).coefs().col(1) -= etaUpdate.patch(p).coefs();
                }
            gsMultiPatch<> geoALEnow;
            for (index_t p = 0; p < 3; ++p)
                geoALEnow.addPatch(nsAssembler.patches().patch(p+3).clone());
            geoALEnow.computeTopology();
            gsElPoissonAssembler<real_t> aleAssembler(geoALEnow,basisALE,bcInfoALE,gALE);
            aleAssembler.options().setReal("LocalStiff",meshStiff);
            aleAssembler.setFixedDofs(0,boundary::south,dispDiff.patch(0).boundary(boundary::north)->coefs().col(0));
            aleAssembler.setFixedDofs(1,boundary::north,dispDiff.patch(0).boundary(boundary::south)->coefs().col(0));
            aleAssembler.setFixedDofs(2,boundary::west,dispDiff.patch(0).boundary(boundary::east)->coefs().col(0));
            aleAssembler.assemble();
    #ifdef GISMO_WITH_PARDISO
            gsSparseSolver<>::PardisoLDLT solverXi(aleAssembler.matrix());
            gsVector<> solVectorXi = solverXi.solve(aleAssembler.rhs());
    #else
            gsSparseSolver<>::SimplicialLDLT solverXi(aleAssembler.matrix());
            gsVector<> solVectorXi = solverXi.solve(aleAssembler.rhs());
    #endif

            // construct ALE update
            aleAssembler.constructSolution(solVectorXi,xiUpdate);

            aleAssembler.setFixedDofs(0,boundary::south,dispDiff.patch(0).boundary(boundary::north)->coefs().col(1));
            aleAssembler.setFixedDofs(1,boundary::north,dispDiff.patch(0).boundary(boundary::south)->coefs().col(1));
            aleAssembler.setFixedDofs(2,boundary::west,dispDiff.patch(0).boundary(boundary::east)->coefs().col(1));
            aleAssembler.assemble();
    #ifdef GISMO_WITH_PARDISO
            gsSparseSolver<>::PardisoLDLT solverEta(aleAssembler.matrix());
            gsVector<> solVectorEta = solverEta.solve(aleAssembler.rhs());
    #else
            gsSparseSolver<>::SimplicialLDLT solverEta(aleAssembler.matrix());
            gsVector<> solVectorEta = solverEta.solve(aleAssembler.rhs());
    #endif

            // construct ALE update
            aleAssembler.constructSolution(solVectorEta,etaUpdate);
            for (index_t p = 0; p < 3; ++p)
            {
                ALE.patch(p).coefs().col(0) += xiUpdate.patch(p).coefs();
                ALE.patch(p).coefs().col(1) += etaUpdate.patch(p).coefs();
                nsAssembler.patches().patch(p+3).coefs().col(0) += xiUpdate.patch(p).coefs();
                nsAssembler.patches().patch(p+3).coefs().col(1) += etaUpdate.patch(p).coefs();
            }

            // velocity ALE
            gsMultiPatch<> velocityALE;
            for (index_t p = 0; p < 3; ++p)
            {
                velocityALE.addPatch(geoALE.patch(p).clone());
                velocityALE.patch(p).coefs().col(0) = xiUpdate.patch(p).coefs()/timeStep;
                velocityALE.patch(p).coefs().col(1) = etaUpdate.patch(p).coefs()/timeStep;
            }
            timeALE += clock.stop();


            // FLOW
            clock.restart();

            nsAssembler.setFixedDofs(3,boundary::south,velocityALE.patch(0).boundary(boundary::south)->coefs());
            nsAssembler.setFixedDofs(4,boundary::north,velocityALE.patch(1).boundary(boundary::north)->coefs());
            nsAssembler.setFixedDofs(5,boundary::west,velocityALE.patch(2).boundary(boundary::west)->coefs());
            if (iter > 0)
                nsTimeSolver.recoverState();
            nsTimeSolver.makeTimeStepFSI2(timeStep,velocityALE,alePatches);
            nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
            timeFlow += clock.stop();

            std::vector<std::pair<index_t, boxSide> > bdrySides;
            bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
            bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
            bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
            bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
            bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
            bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
            gsMatrix<> force = nsAssembler.computeForce(velocity,pressure,bdrySides);
            gsInfo << "Drag: " << force.at(0) << std::endl;
            gsInfo << "Lift: " << force.at(1) << std::endl;

            ++iter;
        }

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,i+1,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,i+1);
        }
        validation(file,nsAssembler,(i+1)*timeStep,velocity,pressure,displacement);
    }

    gsInfo << "Time ALE: " << timeALE << "s, time flow: " << timeFlow << "s, time beam: " << timeBeam << "s.\n";
    //=============================================//
                   // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        gsInfo << "Open \"flappingBeam_FSI2_*.pvd\" in Paraview for visualization.\n";
        collectionFlow.save();
        collectionBeam.save();
        collectionALE.save();
    }
    file.close();

    return 0;
}
