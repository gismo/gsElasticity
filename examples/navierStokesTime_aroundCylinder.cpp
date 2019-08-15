/// This is an example of using the time-dependent Navier-Stokes solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){
    gsInfo << "Testing the time-dependent Navier-Stokes solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t viscosity = 0.001;
    real_t density = 1.;
    real_t maxInflow = 1.5;
    bool subgrid = false;
    bool supg = false;
    index_t numTimeSteps = 100;
    real_t timeSpan = 10.;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the time-dependent Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addInt("s","steps","Number of time steps",numTimeSteps);
    cmd.addReal("d","density","Density of the fluid",density);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,&inflow,0);
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
    }

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisVelocity(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    gsNsAssembler<real_t> stiffAssembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    stiffAssembler.options().setReal("Viscosity",viscosity);
    stiffAssembler.options().setInt("DirichletValues",dirichlet::interpolation);
    stiffAssembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized system with " << stiffAssembler.numDofs() << " dofs.\n";

    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    gsNsTimeIntegrator<real_t> timeSolver(stiffAssembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    timeSolver.options().setInt("Verbosity",newton_verbosity::all);
    // set initial conditions
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
    gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,0,numPlotPoints);

    real_t timeStep = timeSpan/numTimeSteps;
    gsStopwatch clock;
    gsProgressBar bar;

    //=============================================//
              // Transient simulation //
    //=============================================//

    clock.restart();
    gsInfo << "Running the transient simulation...\n";
    for (index_t i = 0; i < numTimeSteps; ++i)
    {
        bar.display(i+1,numTimeSteps);
        timeSolver.makeTimeStep(timeStep);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";
    collection.save();
    gsInfo << "The results are saved to the Paraview file \"NS_aroundCylinder.pvd\"\n";

    return 0;
}
