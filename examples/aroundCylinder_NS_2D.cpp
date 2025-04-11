/// This is the incompressible Navier-Stokes solver benchmark 2D-1 from this project:
/// http://www.featflow.de/en/benchmarks/cfdbenchmarking.html
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsNsAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

void refineBoundaryLayer(gsMultiBasis<> & velocity, gsMultiBasis<> & pressure)
{
    gsMatrix<> boxSouth(2,2);
    boxSouth << 0.,0.,0.,0.2;
    for (index_t p = 0; p < 4; ++p)
    {
        velocity.refine(p,boxSouth);
        pressure.refine(p,boxSouth);
    }
}

int main(int argc, char* argv[]){
    gsInfo << "Benchmark 2D-1: steady-state flow of an incompressible fluid\n";

    //=====================================//
                    // Input //
    //=====================================//

    std::string filename = gsElasticity_DATA"/flow_around_cylinder.xml";
    real_t viscosity = 0.001; // kinematic viscosity
    real_t meanVelocity = 0.2; // inflow velocity
    real_t density = 1.;
    // space discretization
    index_t numUniRef = 1; //was 3
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

    // inflow velocity profile U(y) = U_max*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

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

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
                      // Solving //
    //=============================================//

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setReal("Density",density);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // set assembly type: linear system for Newton's method in the next-solution form
    assembler.options().setInt("Assembly",ns_assembly::newton_next);
    // setting Newton's method
    gsIterative<real_t> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setInt("Solver",linear_solver::LU);
    // set iterative solver mode: each iteration yields an update to the solution
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
        gsWriteParaviewMultiPhysics(fields,"aroundCylinderNS",numPlotPoints,plotMesh);
        gsInfo << "Open \"aroundCylinderNS.pvd\" in Paraview for visualization.\n";
    }

    // computing forces acting on the surface of the solid body
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);

    // evaluating pressure difference at the far front and the far rear points of the cylinder
    gsMatrix<> point(2,1);
    point << 0.5, 0;
    real_t pressureFront = pressure.patch(0).eval(point)(0,0);
    real_t pressureBack = pressure.patch(2).eval(point)(0,0);
    real_t L = 0.1; // characteristic length

    gsInfo << "Drag coefficient: " << 2.*force.at(0)/L/pow(meanVelocity,2) << std::endl;
    gsInfo << "Lift coefficient: " << 2.*force.at(1)/L/pow(meanVelocity,2) << std::endl;
    gsInfo << "Pressure difference: " << pressureFront - pressureBack << std::endl;

    return 0;
}
