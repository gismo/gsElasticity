/// This is the incompressible Stokes solver benchmark based on this project:
/// http://www.featflow.de/en/benchmarks/cfdbenchmarking.html
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
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
    gsInfo << "Benchmark CFD2D0: steady-state flow of an incompressible super-viscous fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numPlotPoints = 10000;
    real_t viscosity = 0.001;
    real_t density = 1.;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    bool plotMesh = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CFD2D0: steady-state flow of an incompressible super-viscous fluid.");
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
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();


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
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
#endif

    // constructing solution as an IGA function
    gsMultiPatch<> velocity, pressure;
    assembler.constructSolution(solVector,velocity,pressure);

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
    gsWriteParaviewMultiPhysics(fields,"aroundCylinder_CFD2D0",numPlotPoints,plotMesh);
    gsInfo << "Open \"aroundCylinder_CFD2D0.pvd\" in Paraview for visualization.\n";

    return 0;
}
