/// This is the incompressible Navier-Stokes solver benchmark CFD1 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
///
/// In this version, the linearization of the INSE is used that yields a new solution (not an update) at each iteration.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsNsAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

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
    gsInfo << "Benchmark CFD1: steady-state flow of an incompressible fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = gsElasticity_DATA"/flappingBeam_flow.xml";
    real_t viscosity = 0.001; // kinematic viscosity
    real_t meanVelocity = 0.2; // inflow velocity
    real_t density = 1.0e3;
    // space discretization
    index_t numUniRef = 3;
    index_t numDegElev = 0;
    index_t numBLRef = 1;
    bool subgridOrTaylorHood = false;
    // output
    index_t numPlotPoints = 10000;
    bool plotMesh = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevations",numDegElev);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addSwitch("e","element","Mixed element: false = subgrid (default), true = Taylor-Hood",subgridOrTaylorHood);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("m","mesh","Plot computational mesh",plotMesh);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisVelocity(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    // additional velocity refinement for stable mixed FEM
    if (!subgridOrTaylorHood) // subgrid
        basisVelocity.uniformRefine();
    else // Taylor-Hood
        basisVelocity.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0); // x-component
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1); // y-component
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(5,boundary::west,condition_type::dirichlet,0,d);
        bcInfo.addCondition(6,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(6,boundary::north,condition_type::dirichlet,0,d);
    }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
                      // Solving //
    //=============================================//

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setReal("Density",density);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    assembler.options().setInt("Assembly",ns_assembly::newton_next);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative<real_t> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setInt("Solver",linear_solver::LU);
    solver.options().setInt("IterType",iteration_type::next);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    solver.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                      // Output //
    //=============================================//

    // solution as two isogeometric fields
    gsMultiPatch<> velocity, pressure;
    assembler.constructSolution(solver.solution(),solver.allFixedDofs(),velocity,pressure);

    if (numPlotPoints > 0) // visualization
    {
        // constructing isogeometric field (geometry + solution)
        gsField<> velocityField(assembler.patches(),velocity);
        gsField<> pressureField(assembler.patches(),pressure);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Velocity"] = &velocityField;
        fields["Pressure"] = &pressureField;
        gsWriteParaviewMultiPhysics(fields,"flappingBeam_CFD1",numPlotPoints,plotMesh);
        gsInfo << "Open \"flappingBeam_CFD1.pvd\" in Paraview for visualization.\n";
    }

    // computing forces acting on the surface of the solid body
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);
    gsInfo << "Drag: " << force.at(0) << std::endl;
    gsInfo << "Lift: " << force.at(1) << std::endl;

    return 0;
}
