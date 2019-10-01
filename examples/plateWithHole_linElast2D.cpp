/// This is the 2D linear elasticity benchmark "Infinite plate with circular hole"
/// as described in Hughes, et.al. 2005
/// "Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement"
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is the 2D linear elasticity benchmark: infinite plate with circular hole.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/plateWithHole.xml";
    index_t numUniRef = 5; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 10000;
    bool plotMesh = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is the 2D linear elasticity benchmark: infinite plate with circular hole.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("m","mesh","Plot computational mesh",plotMesh);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating basis
    gsMultiBasis<> basis(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basis.degreeElevate();
        basis.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // boundary load neumann BC
    gsFunctionExpr<> traction("-10*(x==-4)","10*(y==4)",2);
    //gsConstantFunction<> traction(10.,0.,2);

    // material parameters
    real_t youngsModulus = 1e5;
    real_t poissonsRatio = 0.3;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::neumann,&traction);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
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
    gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
#endif

    // constructing solution as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,solution);
    gsInfo << "Checking bijectivity...\n";
    if (assembler.checkSolution(solution) != -1)
        gsInfo << "Computed displacement field is not valid (J < 0)!\n";

    // constructing an IGA field (geometry + solution)
    gsField<> solutionField(assembler.patches(),solution);

    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_type::all_2D);
    gsField<> stressField(assembler.patches(),stresses,true);

    //=============================================//
                  // Output //
    //=============================================//

    if (numPlotPoints > 0)
    {
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Stresses"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"plateWithHole",numPlotPoints,plotMesh);
        gsInfo << "Open \"plateWithHole.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
