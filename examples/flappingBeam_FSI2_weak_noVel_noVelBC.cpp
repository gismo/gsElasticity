/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
///
/// weak coupling with new geo
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

void refineBoundaryLayer(gsMultiBasis<> & velocity, gsMultiBasis<> & pressure)
{
    gsMatrix<> boxEast(2,2);
    boxEast << 0.8,1.,0.,0.;
    velocity.refine(0,boxEast);
    pressure.refine(0,boxEast);

    gsMatrix<> boxSouth(2,2);
    boxSouth << 0.,0.,0.,0.2;
    velocity.refine(1,boxSouth);
    velocity.refine(3,boxSouth);
    pressure.refine(1,boxSouth);
    pressure.refine(3,boxSouth);

    gsMatrix<> boxNorth(2,2);
    boxNorth << 0.,0.,0.8,1.;
    velocity.refine(2,boxNorth);
    velocity.refine(4,boxNorth);
    pressure.refine(2,boxNorth);
    pressure.refine(4,boxNorth);

    gsMatrix<> boxWest(2,2);
    boxWest << 0.,0.2,0.,0.;
    velocity.refine(5,boxWest);
    pressure.refine(5,boxWest);
}

void formResidual(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow,
                  gsMatrix<> & residual)
{
    index_t numBdry = interfaceOld.size();
    short_t dim = interfaceOld[0].cols();
    index_t numDoFs = 0;
    for (index_t i = 0; i < numBdry; ++i)
        numDoFs += dim*interfaceOld[i].rows();
    residual.setZero(numDoFs,1);

    index_t index = 0;
    for (index_t i = 0; i < numBdry; ++i)
    {
        index_t size = interfaceOld[i].rows();
        for (index_t d = 0; d < dim; ++d)
        {
            residual.middleRows(index,size) = interfaceNow[i].col(d) - interfaceOld[i].col(d);
            index += size;
        }
    }
}

void aitkenRelaxation(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow,
                      std::vector<gsMatrix<> > & interfaceNew, real_t & omega)
{
    gsMatrix<> residualOld,residualNew;
    formResidual(interfaceOld,interfaceNow,residualOld);
    formResidual(interfaceNow,interfaceNew,residualNew);

    omega = -1.*omega*(residualOld.transpose()*(residualNew-residualOld))(0,0) /
            pow((residualNew-residualOld).norm(),2);

    for (index_t i = 0; i < interfaceOld.size(); ++i)
    {
        interfaceOld[i] = interfaceNow[i];
        interfaceNow[i] += omega*(interfaceNew[i]-interfaceNow[i]);
    }
}

real_t computeResidual(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow)
{
    gsMatrix<> residual;
    formResidual(interfaceOld,interfaceNow,residual);
    return residual.norm();
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the steady fluid-structure interaction solver in 2D.\n";

    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameFlowPart = ELAST_DATA_DIR"/flappingBeam_flowPart2.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRefFlow = 3; // number of h-refinements for the fluid
    index_t numKRefFlow = 0; // number of k-refinements for the fluid
    index_t numBLRef = 0; // number of additional boundary layer refinements for the fluid
    index_t numUniRefBeam = 3; // number of h-refinements for the beam and the ALE mapping
    index_t numKRefBeam = 0; // number of k-refinements for the beam and the ALE mapping
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;//4.4e6;//1.4e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    bool subgrid = true;
    bool supg = false;
    real_t densityFluid = 1.0e3;
    real_t densitySolid = 1.0e4;
    real_t absTol = 1e-10;
    real_t relTol = 1e-10;
    real_t timeStep = 0.001;
    real_t timeSpan = 0.001;
    real_t warmUpTimeSpan = 3.;
    real_t warmUpTimeStep = 0.1;

    real_t meshPR = 0.4;
    real_t stiff = 3.;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the steady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRefFlow);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("b","beamrefine","Number of uniform refinement applications for the beam and ALE",numUniRefBeam);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Pois ration for mesh",meshPR);
    cmd.addReal("j","stiff","",stiff);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoPart; // this is a part of the flow geometry; we deform only this part to save memory and time
    gsReadFile<>(filenameFlowPart, geoPart);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);
    geoPart.computeTopology();

    // creating bases
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisPressure(geoFlow);
    for (index_t i = 0; i < numKRefFlow; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRefFlow; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);

    // additional velocity refinement for stable mixed FEM;
    // displacement is also refined to keep bases match at the fsi interface
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numKRefBeam; ++i)
    {
        basisDisplacement.degreeElevate();
        geoPart.degreeElevate();
        geoFlow.degreeElevate();
        basisDisplacement.uniformRefine();
        geoPart.uniformRefine();
        geoFlow.uniformRefine();
    }
    for (index_t i = 0; i < numUniRefFlow; ++i)
    {
        basisDisplacement.uniformRefine();
        geoPart.uniformRefine();
        geoFlow.uniformRefine();
    }
    basisDisplacement.uniformRefine();
    geoPart.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisALE(geoPart);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
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
        /*bcInfoFlow.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(5,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(7,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(9,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(9,boundary::north,condition_type::dirichlet,0,d);*/
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
    gsFsiLoad<real_t> fSouth(geoPart,ALE,4,boundary::north,
                             velocity,pressure,4,viscosity,densityFluid);
    gsFsiLoad<real_t> fEast(geoPart,ALE,5,boundary::west,
                            velocity,pressure,5,viscosity,densityFluid);
    gsFsiLoad<real_t> fNorth(geoPart,ALE,3,boundary::south,
                             velocity,pressure,3,viscosity,densityFluid);
    bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&fSouth);
    bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&fEast);
    bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&fNorth);
    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoALE;
    for (auto it = geoPart.bBegin(); it != geoPart.bEnd(); ++it)
        for (index_t d = 0; d < 2; ++d)
            bcInfoALE.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

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
    nsTimeSolver.options().setInt("Scheme",time_integration_NS::theta_scheme_linear);
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
    gsElasticityAssembler<real_t> aleAssembler(geoPart,basisALE,bcInfoALE,g);
    aleAssembler.options().setReal("PoissonsRatio",meshPR);
    aleAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    aleAssembler.options().setReal("LocalStiff",stiff);
    gsInfo << "Initialized elasticity system for ALE with " << aleAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velocity);
    gsField<> pressureField(nsAssembler.patches(),pressure);
    gsField<> displacementField(geoBeam,displacement);
    gsField<> aleField(geoPart,ALE);

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

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());
    nsTimeSolver.initialize();

    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    for (index_t i = 0; i < index_t(warmUpTimeSpan/warmUpTimeStep); ++i)
    {
        bar.display(i+1,index_t(warmUpTimeSpan/warmUpTimeStep));
        if ((i+1)*warmUpTimeStep < 2./3.*warmUpTimeSpan)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*warmUpTimeStep*(i+1)/2./warmUpTimeSpan*3.))/2);
        else
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs);
        nsTimeSolver.makeTimeStep(warmUpTimeStep);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
             // Initial conditions //
    //=============================================//

    // plotting initial condition
    nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
    elAssembler.constructSolution(gsMatrix<>::Zero(elAssembler.numDofs(),1),displacement);
    aleAssembler.constructSolution(gsMatrix<>::Zero(aleAssembler.numDofs(),1),ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_test_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_test_beam",collectionBeam,0,numPlotPoints);
        //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_test_ALE",collectionALE,0,numPlotPoints);
        //plotDeformation(geoPart,ALE,"flappingBeam_FSI2_test_ALE",collectionALE,0);
    }

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

    std::vector<gsMatrix<> > interfaceNow, interfaceAitken, interfaceOld, interfaceNew;
    interfaceNow.push_back(ALE.patch(3).boundary(boundary::south)->coefs());
    interfaceNow.push_back(ALE.patch(4).boundary(boundary::north)->coefs());
    interfaceNow.push_back(ALE.patch(5).boundary(boundary::west)->coefs());
    interfaceNew = interfaceNow;


    std::vector<std::pair<index_t,index_t> > alePatches;
    alePatches.push_back(std::pair<index_t,index_t>(0,0));
    alePatches.push_back(std::pair<index_t,index_t>(1,1));
    alePatches.push_back(std::pair<index_t,index_t>(2,2));
    alePatches.push_back(std::pair<index_t,index_t>(3,3));
    alePatches.push_back(std::pair<index_t,index_t>(4,4));
    alePatches.push_back(std::pair<index_t,index_t>(5,5));
    //alePatches.push_back(std::pair<index_t,index_t>(6,6));
    //alePatches.push_back(std::pair<index_t,index_t>(7,7));
    //alePatches.push_back(std::pair<index_t,index_t>(8,8));
    //alePatches.push_back(std::pair<index_t,index_t>(9,9));

    gsNewton<real_t> aleNewton(aleAssembler,gsMatrix<>::Zero(aleAssembler.numDofs(),1),aleAssembler.allFixedDofs());
    aleNewton.options().setInt("Verbosity",newton_verbosity::none);
    aleNewton.options().setInt("MaxIters",1);
    aleNewton.options().setInt("Solver",linear_solver::LDLT);

    elTimeSolver.initialize();

    gsMultiPatch<> velocityALE;
    aleAssembler.constructSolution(gsMatrix<>::Zero(aleAssembler.numDofs(),1),aleAssembler.allFixedDofs(),velocityALE);

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    clock.restart();
    gsInfo << "Running the coupled simulation...\n";
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        gsInfo << "=========================================================TIME STEP " << i+1 << "/" << index_t(timeSpan/timeStep) << std::endl;

        if (aleAssembler.checkSolution(ALE) != -1)
            break;


        // 1. FLOW
        // set ALE velocity on the interface
        //nsAssembler.setFixedDofs(3,boundary::south,velocityALE.patch(3).boundary(boundary::south)->coefs());
        //nsAssembler.setFixedDofs(4,boundary::north,velocityALE.patch(4).boundary(boundary::north)->coefs());
        //nsAssembler.setFixedDofs(6,boundary::west,velocityALE.patch(6).boundary(boundary::west)->coefs());
        // give all info to the flow time solver to make a time step
        //nsTimeSolver.makeTimeStep(timeStep);

        nsTimeSolver.makeTimeStepFSI2(timeStep,velocityALE,alePatches);
        // construct velocity and pressure
        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),nsTimeSolver.allFixedDofs(),velocity,pressure);
        //velocity.patch(3).coefs() -= velocityALE.patch(0).coefs();
        //velocity.patch(4).coefs() -= velocityALE.patch(1).coefs();
        //velocity.patch(5).coefs() -= velocityALE.patch(2).coefs();
        gsMatrix<> force = nsAssembler.computeForce(velocity,pressure,bdrySides);
        gsInfo << "Drag: " << force.at(0) << std::endl;
        gsInfo << "Lift: " << force.at(1) << std::endl;

        // 2. BEAM
        elTimeSolver.makeTimeStep(timeStep);
        elAssembler.constructSolution(elTimeSolver.displacementVector(),displacement);

        interfaceNew[0] = displacement.patch(0).boundary(boundary::north)->coefs();
        interfaceNew[1] = displacement.patch(0).boundary(boundary::south)->coefs();
        interfaceNew[2] = displacement.patch(0).boundary(boundary::east)->coefs();

        // 4. ALE
        // compute ALE velocity
        aleAssembler.constructSolution(aleNewton.solution(),aleNewton.allFixedDofs(),velocityALE);

        // reverse previous deformation of the flow domain (we will overwrite ALE)
        for (index_t p = 0; p < 6; ++p)
            nsAssembler.patches().patch(p).coefs() -= ALE.patch(p).coefs();
        // set DBC for ALE update: new Aitken interface minus interface at the start of the time step
        aleAssembler.setFixedDofs(3,boundary::south,interfaceNew[0]-interfaceNow[0]);
        aleAssembler.setFixedDofs(4,boundary::north,interfaceNew[1]-interfaceNow[1]);
        aleAssembler.setFixedDofs(5,boundary::west,interfaceNew[2]-interfaceNow[2]);
        // compute ALE update; ALE at the start of the time step is used as a current configuration
        aleNewton.reset();
        aleNewton.solve();
        // construct new ALE
        aleAssembler.constructSolution(aleNewton.solution(),aleNewton.allFixedDofs(),ALE);
        // apply new deformation to the flow domain
        for (index_t p = 0; p < 6; ++p)
            nsAssembler.patches().patch(p).coefs() += ALE.patch(p).coefs();
        //for (index_t p = 0; p < 3; ++p)
        //    velocityALE.patch(p).coefs() = (ALE.patch(p).coefs() - velocityALE.patch(p).coefs()) / timeStep;

        interfaceNow = interfaceNew;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_test_flow",collectionFlow,i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_test_beam",collectionBeam,i+1,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_test_ALE",collectionALE,i+1,numPlotPoints);
            //plotDeformation(geoPart,ALE,"flappingBeam_FSI2_test_ALE",collectionALE,i+1);
        }

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
        //collectionALE.save();
    }

    return 0;
}
