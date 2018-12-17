/// This is an example of using the linear elasticity solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the linear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/lshape.xml";
    int numUniRef = 3; // number of h-refinements
    int numDegElevate = 1; // number of p-refinements
    int numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","prefine","Number of degree elevation application",numDegElevate);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // boundary displacement in y-direction, dirichlet BC
    gsFunctionExpr<> bdry_disp("0.1",2);

    // material parameters
    real_t youngsModulus = 200.;//74e9;
    real_t poissonsRatio = 0.33;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1);
    bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,0,0);
    bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,&bdry_disp,1);

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basis(geometry);
    for (int i = 0; i < numDegElevate; ++i)
        basis.degreeElevate();
    for (int i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    gsInfo<<"Assembling...\n";
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
    gsMatrix<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with LU solver.\n";

    // constructing solution as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,solution,gsVector<index_t>::vec(0,1));
    gsInfo << "Constructed solution.\n";
    // constructing an IGA field (geometry + solution)
    gsField<> solutionField(assembler.patches(),solution);

    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_type::von_mises);
    gsField<> stressField(assembler.patches(),stresses,true);

    //=============================================//
                  // Output //
    //=============================================//

    gsInfo << "Plotting the output to Paraview...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Deformation"] = &solutionField;
    fields["Stresses"] = &stressField;
    gsWriteParaviewMultiPhysics(fields,"lshape",numPlotPoints);
    gsInfo << "Finished.\n";
    gsInfo << "Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
