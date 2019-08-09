/// This is an example of using the Navier-Stokes solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){
    gsInfo << "Testing the Navier-Stokes solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t viscosity = 0.001;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    index_t iters = 50;
    bool supg = true;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("i","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addInt("j","iters","Max number of Newton's iterations",iters);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
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
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    assembler.options().setSwitch("SUPG",supg);
    assembler.options().setSwitch("Iteration",false);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving Oseen//
    //=============================================//

    gsMultiPatch<> velocity, pressure;

    gsField<> velocityField(assembler.patches(),velocity);
    gsField<> pressureField(assembler.patches(),pressure);

    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;

    gsParaviewCollection collection("NS_around_cylinder");


    index_t numIter = iters;
    assembler.assemble();
    gsSparseSolver<>::LU solver(assembler.matrix());
    gsMatrix<> solVector = solver.solve(assembler.rhs());
    assembler.constructSolution(solVector,velocity,pressure);
    gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_around_cylinder",collection,0,numPlotPoints);

    gsMatrix<> tempSolVector;
    for (index_t i = 0; i < 6; ++i)
    {
        assembler.assemble(solVector);
        gsSparseSolver<>::LU solver(assembler.matrix());
        tempSolVector = solver.solve(assembler.rhs());
        gsInfo << "It " << i+1 << " abs " << (tempSolVector-solVector).norm() << std::endl;
        solVector = tempSolVector;
        assembler.constructSolution(solVector,velocity,pressure);
        gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_around_cylinder",collection,i+1,numPlotPoints);

    }

    //=============================================//
                  // Solving Newton//
    //=============================================//

    assembler.options().setSwitch("Iteration",true);
    // setting Newton's method
    gsNewton<real_t> newton(assembler,solVector);
    newton.options().setInt("Verbosity",newton_verbosity::all);
    newton.options().setInt("MaxIters",iters);
    newton.options().setInt("Solver",linear_solver::LU);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    // constructing solution as an IGA function
    //gsMultiPatch<> velocity, pressure;
    assembler.constructSolution(newton.solution(),velocity,pressure);
    gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_around_cylinder",collection,7,numPlotPoints);

    // constructing an IGA field (geometry + solution)

    //=============================================//
                  // Output //
    //=============================================//

    gsInfo << "Plotting the output to the Paraview file \"NS_around_cylinder.pvd\"...\n";
    // creating a container to plot all fields to one Paraview file
    //gsWriteParaviewMultiPhysics(fields,"NS_around_cylinder",numPlotPoints);

    //=============================================//
                  // Validation //
    //=============================================//

    gsMatrix<> point(2,1);
    point << 0.5, 0;
    gsInfo << "Pressure difference: " << pressure.patch(0).eval(point)(0,0) -
                                         pressure.patch(2).eval(point)(0,0) << "Pa\n";

    collection.save();

    return 0;

}
