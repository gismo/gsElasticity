/// This is an example of using the linear elasticity solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the linear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = "lshape.xml";
    int numUniRef = 3; // number of h-refinements
    int numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> f(0.,0.,2);
    /*
    // boundary displacement in y-direction, dirichlet BC
    gsConstantFunction<> g(0.1,2);
    */
    gsConstantFunction<> g(0.0, 0.1,2);

    // material parameters
    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;
    real_t density = 2.8e3; // doesn't really matter for a stationary problem

    // boundary conditions
    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(2,boundary::east,condition_type::dirichlet,0); // third number is a component (coordinate) number
    bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,&g);

/*
    // Dirichlet BC are imposed separately for every component (coordinate)
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0);  // first number is a patch number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1);  // second number is actually a NULL pointer, meaning 0 function
    bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,0,0); // third number is a component (coordinate) number
    bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,&g,1);
*/

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
    gsElasticityAssembler<real_t> assembler(geometry,basis,youngsModulus,poissonsRatio,density,bcInfo,f);
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
    assembler.constructSolution(solVector,solution);
    // constructing an IGA field (geometry + solution)
    gsField<> solutionField(assembler.patches(),solution);

    // constructing stresses
    gsMultiFunction<real_t> vonMisesStresses;
    assembler.constructStresses(solVector,vonMisesStresses,stress_type::von_mises);
    gsField<> vonMisesStressField(assembler.patches(),vonMisesStresses,true);


    //=============================================//
                  // Output //
    //=============================================//

    gsInfo << "Internal energy: " << solVector.transpose() * assembler.matrix() * solVector << "\n";
    gsInfo << "External energy: " << assembler.rhs().transpose() * solVector << "\n";

    gsInfo << "Plotting the output to Paraview...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Deformation"] = &solutionField;
    fields["vonMises"] = &vonMisesStressField;
    gsWriteParaviewMultiPhysics(fields,"lshape",numPlotPoints);
    gsInfo << "Finished.\n";
    gsInfo << "Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
