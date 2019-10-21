/// This is the incompressible Navier-Stokes solver benchmark 2D-1 from this project:
/// http://www.featflow.de/en/benchmarks/cfdbenchmarking.html
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

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
    gsInfo << "Benchmark CFD2D1: steady-state flow of an incompressible fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numPlotPoints = 10000;
    bool plotMesh = false;
    real_t viscosity = 0.001;
    real_t density = 1.;
    real_t maxInflow = 0.3;
    bool subgrid = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("m","mesh","Plot computational mesh",plotMesh);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
                  // Assembly //
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
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    // additional velocity refinement for stable mixed FEM
    if (subgrid)
    {
        gsInfo << "Using subgrid element.\n";
        basisVelocity.uniformRefine();
    }
    else
    {
        gsInfo << "Using Taylor-Hood element.\n";
        basisVelocity.degreeElevate();
    }

    // inflow velocity profile U(y) = U_max*y*(H-y)/(H/2)^2; channel height H = 0.41
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

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setReal("Density",density);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    assembler.options().setInt("Iteration",iteration_type::newton);
    // setting Newton's method
    gsNewton<real_t> newton(assembler);
    newton.options().setInt("Verbosity",newton_verbosity::all);
    newton.options().setInt("Solver",linear_solver::LU);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    // solution as two isogeometric fields
    gsMultiPatch<> velocity, pressure;
    assembler.constructSolution(newton.solution(),newton.allFixedDofs(),velocity,pressure);

    //=============================================//
                  // Validation //
    //=============================================//

    // computing forces acting on the surface of the solid body
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);

    real_t L = 0.1; // characteristic lenght
    gsInfo << "Drag: " << force.at(0) << std::endl;
    gsInfo << "Drag coefficient: " << 2.*force.at(0)/L/pow(maxInflow*2./3.,2) << std::endl;
    gsInfo << "Lift: " << force.at(1) << std::endl;
    gsInfo << "Lift coefficient: " << 2.*force.at(1)/L/pow(maxInflow*2./3.,2) << std::endl;
    // evaluating pressure difference at the far front and the far rear points of the cylinder
    gsMatrix<> point(2,1);
    point << 0.5, 0;                                                                   // this info is hard-cored in the geometry
    gsInfo << "Pressure difference: " << pressure.patch(0).eval(point)(0,0) -          // far front point
                                         pressure.patch(2).eval(point)(0,0) << "Pa\n"; // far rear point

    //=============================================//
                  // Visualization //
    //=============================================//

    // constructuin isogeometric field (geometry + solution)
    gsField<> velocityField(assembler.patches(),velocity);
    gsField<> pressureField(assembler.patches(),pressure);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;
    gsWriteParaviewMultiPhysics(fields,"aroundCylinder_CFD2D1",numPlotPoints,plotMesh);
    gsInfo << "Open \"aroundCylinder_CFD2D1.pvd\" in Paraview for visualization.\n";

    return 0;
}
