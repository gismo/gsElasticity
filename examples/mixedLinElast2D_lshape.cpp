/// This is an example of using the mixed linear elasticity solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsMixedElasticityAssembler.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the linear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/lshape.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numDegElevate = 1; // number of p-refinements
    index_t numPlotPoints = 10000;

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
    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElevate; ++i)
    {
        basisDisplacement.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    basisDisplacement.degreeElevate();

    // creating assembler
    gsMixedElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
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
    gsSparseSolver<>::LU solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with LU solver in " << clock.stop() <<"s.\n";

    // constructing solution as an IGA function
    gsMultiPatch<> displacement, pressure;
    assembler.constructSolution(solVector,displacement,pressure);

    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(),displacement);
    gsField<> pressureField(assembler.patches(),pressure);

    //=============================================//
                  // Output //
    //=============================================//

    gsInfo << "Plotting the output to the Paraview file \"lshape.pvd\"...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &displacementField;
    fields["Pressure"] = &pressureField;
    gsWriteParaviewMultiPhysics(fields,"lshape",numPlotPoints);
    gsInfo << "Done. Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
