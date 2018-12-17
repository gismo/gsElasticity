/// This is an example of using the linear elasticity solver on a 3D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the linear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"terrific.xml";
    int numUniRef = 0; // number of h-refinements
    int numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> f(0.,0.,0.,3);
    // surface load, neumann BC
    gsConstantFunction<> g(20e6, -14e6, 0,3);

    // material parameters
    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BC are imposed separately for every component (coordinate)
    for (int d = 0; d < 3; d++)
    {
        bcInfo.addCondition(0,boundary::back,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::back,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
    }
    // Neumann BC are imposed as one function
    bcInfo.addCondition(13,boundary::front,condition_type::neumann,&g);
    bcInfo.addCondition(14,boundary::north,condition_type::neumann,&g);

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basis(geometry);
    for (int i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,f);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);
    gsInfo<<"Assembling...\n";
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    gsSparseSolver<>::LU solver(assembler.matrix());
    gsMatrix<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with LU solver.\n";

    // constructing solution as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,solution,gsVector<index_t>::vec(0,1,2));
    // constructing an IGA field (geometry + solution)
    gsField<> solutionField(assembler.patches(),solution);

    // constructing stresses
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
    gsWriteParaviewMultiPhysics(fields,"terrific",numPlotPoints);
    gsInfo << "Finished.\n";
    gsInfo << "Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
