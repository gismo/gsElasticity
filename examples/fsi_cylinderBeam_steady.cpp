/// This is an example of solving the steady 2D Fluid-Structure Interaction problem
/// using the incompressible Navier-Stokes and nonlinear elasticity solvers
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

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

int main(int argc, char* argv[]){
    gsInfo << "Testing the steady fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filenameFlow = ELAST_DATA_DIR"/fsi_flow_around_cylinder.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/fsi_beam_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements
    index_t numPlotPoints = 10000;
    real_t youngModulus = 0.5e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    bool supg = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the time-dependent Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("y","young","Young modulus of the beam materail",youngModulus);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

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

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

    // creating bases
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisPressure(geoFlow);
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisDisplacement.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    // additional velocity refinement for stable mixed FEM;
    // displacement is also refined to keep bases match at the fsi interface
    if (subgrid)
    {
        basisVelocity.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    else
    {
        basisVelocity.degreeElevate();
        basisDisplacement.degreeElevate();
    }

    /// todo //////////////////////////


    gsNsAssembler<real_t> stiffAssembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    stiffAssembler.options().setReal("Viscosity",viscosity);
    stiffAssembler.options().setInt("DirichletValues",dirichlet::interpolation);
    stiffAssembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized system with " << stiffAssembler.numDofs() << " dofs.\n";

    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    gsNsTimeIntegrator<real_t> timeSolver(stiffAssembler,massAssembler);
    timeSolver.setInitialSolution(gsMatrix<>::Zero(stiffAssembler.numDofs(),1));
    timeSolver.initialize();

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // containers for solution as IGA functions
    gsMultiPatch<> velocity, pressure;
    // isogeometric fields (geometry + solution)
    gsField<> velocityField(geometry,velocity);
    gsField<> pressureField(geometry,pressure);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;
    // paraview collection of time steps
    gsParaviewCollection collection("NS_aroundCylinder");
    // plotting initial condition
    stiffAssembler.constructSolution(gsMatrix<>::Zero(stiffAssembler.numDofs(),1),velocity,pressure);
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,0,numPlotPoints);
    // output file for solver validation
    std::ofstream ofs;
    ofs.open ("validation.txt");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
          // Solving for initial conditions //
    //=============================================//


    timeSolver.options().setInt("Scheme",time_integration_NS::implicit_euler_oseen);
    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    for (index_t i = 0; i < 30; ++i)
    {
        bar.display(i+1,30);
        timeSolver.makeTimeStep(0.1);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                  // Main simulation //
    //=============================================//


    timeSolver.options().setInt("Scheme",time_integration_NS::implicit_euler_oseen);
    clock.restart();
    gsInfo << "Running the simulation with a fine time step...\n";
    for (index_t i = 0; i < timeSpan/timeStep; ++i)
    {
        bar.display(i+1,timeSpan/timeStep);
        timeSolver.makeTimeStep(timeStep);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,i+11,numPlotPoints);
        if (validate)
            validation(ofs,stiffAssembler,maxInflow,velocity,pressure);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "The results are saved to the Paraview file \"NS_aroundCylinder.pvd\"\n";
    }
    ofs.close();

    return 0;
}
