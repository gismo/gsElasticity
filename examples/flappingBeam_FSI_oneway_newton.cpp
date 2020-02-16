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
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides,true);
    // print time-drag-lift-pressure-YdispA
    gsMatrix<> A(2,1);
    A << 0.,0.5;
    A = pressure.patch(5).eval(A);
    gsMatrix<> B(2,1);
    B << 1.,0.5;
    B = pressure.patch(0).eval(B);
    gsMatrix<> C(2,1);
    C << 1.,0.5;
    C = displacement.patch(0).eval(C);
    gsMatrix<> D(2,1);
    D << 0.,0.5;
    D = velocity.patch(5).eval(D);

    ofs << time << " " << force(0,0) + force(0,1)  << " " << force(1,0) + force(1,1) << " "<<B.at(0)-A.at(0)
        << " "<< C.at(1) << " " << D.at(1) << " " << force(1,0) << " " << force(1,1) << std::endl;
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the one-way time-dependent fluid-structure interaction solver in 2D.\n";

    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameFlowPart = ELAST_DATA_DIR"/flappingBeam_flowPart.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    bool subgrid = false;
    real_t densityFluid = 1.0e3;
    real_t densitySolid = 1.0e3;
    real_t timeStep = 0.01;
    real_t timeSpan = 1.;
    real_t warmUpTimeSpan = 2.;
    real_t warmUpTimeStep = 0.1;
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.; // local stiffening for ALE
    real_t loading = 2.;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the one-way fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("j","stiff","Local stiffening for ALE",meshStiff);
    cmd.addReal("l","load","Gravity loading acting on the beam",loading);
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
    gsConstantFunction<> gFlow(0.,0.,2);
    gsConstantFunction<> gBeam(0.,loading*densitySolid,2);
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
    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoALE;
    for (auto it = geoALE.bBegin(); it != geoALE.bEnd(); ++it)
        for (index_t d = 0; d < 2; ++d)
            bcInfoALE.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,gFlow);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,gFlow);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler);
    nsTimeSolver.options().setInt("Scheme",time_integration::implicit_linear);
    nsTimeSolver.options().setInt("Verbosity",solver_verbosity::all);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gFlow);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: flow mesh
    gsElasticityAssembler<real_t> aleAssembler(geoALE,basisALE,bcInfoALE,gFlow);
    aleAssembler.options().setReal("PoissonsRatio",meshPR);
    aleAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    aleAssembler.options().setReal("LocalStiff",meshStiff);
    gsIterative<real_t> aleNewton(aleAssembler,gsMatrix<>::Zero(aleAssembler.numDofs(),1),aleAssembler.allFixedDofs());
    aleNewton.options().setInt("Verbosity",solver_verbosity::none);
    aleNewton.options().setInt("MaxIters",1);
    aleNewton.options().setInt("Solver",linear_solver::LDLT);
    gsInfo << "Initialized elasticity system for ALE with " << aleAssembler.numDofs() << " dofs.\n";

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
                   // Warming up //
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
    aleAssembler.constructSolution(gsMatrix<>::Zero(aleAssembler.numDofs(),1),ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,0,numPlotPoints);
        //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,0,numPlotPoints);
        plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,0);
    }
    validation(file,nsAssembler,0.,velocity,pressure,displacement);

    clock.restart();
    index_t numWUS = index_t(warmUpTimeSpan/warmUpTimeStep);
    gsInfo << "Flow warm-up...\n";
    for (index_t i = 0; i < numWUS; ++i)
    {
        bar.display(i+1,numWUS);
        nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*warmUpTimeStep*(i+1)/warmUpTimeSpan))/2);
        nsTimeSolver.makeTimeStep(warmUpTimeStep);

        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,i+1,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,i+1);
        }
        validation(file,nsAssembler,(i+1)*warmUpTimeStep,velocity,pressure,displacement);
    }
    index_t numCDS = 3;
    gsInfo << "Flow calm-down...\n";
    for (index_t i = 0; i < numCDS; ++i)
    {
        bar.display(i+1,numCDS);
        nsTimeSolver.makeTimeStep(timeStep);

        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,numWUS+i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,numWUS+i+1,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,numWUS+i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,numWUS+i+1);
        }
        validation(file,nsAssembler,warmUpTimeSpan + (i+1)*timeStep,velocity,pressure,displacement);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    std::vector<std::pair<index_t,index_t> > alePatches;
    alePatches.push_back(std::pair<index_t,index_t>(3,0));
    alePatches.push_back(std::pair<index_t,index_t>(4,1));
    alePatches.push_back(std::pair<index_t,index_t>(5,2));

    elTimeSolver.initialize();
    clock.restart();
    gsInfo << "FSI...\n";
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));

        if (aleAssembler.checkSolution(ALE) != -1)
            break;

        // BEAM
        gsMultiPatch<> dispDiff;
        elAssembler.constructSolution(elTimeSolver.displacementVector(),dispDiff);
        elTimeSolver.makeTimeStep(timeStep);
        elAssembler.constructSolution(elTimeSolver.displacementVector(),displacement);
        dispDiff.patch(0).coefs() = displacement.patch(0).coefs() - dispDiff.patch(0).coefs();

        // ALE
        aleAssembler.setFixedDofs(0,boundary::south,dispDiff.patch(0).boundary(boundary::north)->coefs());
        aleAssembler.setFixedDofs(1,boundary::north,dispDiff.patch(0).boundary(boundary::south)->coefs());
        aleAssembler.setFixedDofs(2,boundary::west,dispDiff.patch(0).boundary(boundary::east)->coefs());
        aleNewton.reset();
        aleNewton.solve();
        // velocity ALE
        gsMultiPatch<> velocityALE;
        aleAssembler.constructSolution(aleNewton.solution(),aleNewton.allFixedDofs(),velocityALE);
        for (index_t p = 0; p < velocityALE.nPatches(); ++p)
        {
            velocityALE.patch(p).coefs() -= ALE.patch(p).coefs();
            velocityALE.patch(p).coefs() /= timeStep;
        }
        // undo last ALE
        nsAssembler.patches().patch(3).coefs() -= ALE.patch(0).coefs();
        nsAssembler.patches().patch(4).coefs() -= ALE.patch(1).coefs();
        nsAssembler.patches().patch(5).coefs() -= ALE.patch(2).coefs();
        // construct ALE
        aleAssembler.constructSolution(aleNewton.solution(),aleNewton.allFixedDofs(),ALE);
        // apply new ALE
        nsAssembler.patches().patch(3).coefs() += ALE.patch(0).coefs();
        nsAssembler.patches().patch(4).coefs() += ALE.patch(1).coefs();
        nsAssembler.patches().patch(5).coefs() += ALE.patch(2).coefs();


        // FLOW
        nsAssembler.setFixedDofs(3,boundary::south,velocityALE.patch(0).boundary(boundary::south)->coefs());
        nsAssembler.setFixedDofs(4,boundary::north,velocityALE.patch(1).boundary(boundary::north)->coefs());
        nsAssembler.setFixedDofs(5,boundary::west,velocityALE.patch(2).boundary(boundary::west)->coefs());
        nsTimeSolver.makeTimeStepFSI3(timeStep,velocityALE,alePatches);
        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,numWUS+numCDS+i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,numWUS+numCDS+i+1,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,numWUS+numCDS+i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_FSI2_ALE",collectionALE,numWUS+numCDS+i+1);
        }
        validation(file,nsAssembler,warmUpTimeSpan+numCDS*timeStep+(i+1)*timeStep,velocity,pressure,displacement);
    }

    gsInfo << "Complete in " << clock.stop() << "s.\n";
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
