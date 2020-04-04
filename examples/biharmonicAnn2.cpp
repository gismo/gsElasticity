/// This is an example of solving the biharmonic equation using the Poisson's solver.
///
/// Authors: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElPoissonAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the biharmonic solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    index_t numUniRef = 4;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;
    bool subgridOrTaylorHood = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("e","element","Mixed element: false = subgrid (default), true = Taylor-Hood",subgridOrTaylorHood);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );// creating basis
    gsMultiBasis<> basis(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basis.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsFunctionExpr<> source  ("-64*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-4*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)/4",2);
    // surface load, neumann BC

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,&laplace);
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,&laplace);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,&laplace);
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,&laplace);
    gsBoundaryConditions<> bcInfo2;
    bcInfo2.addCondition(0,boundary::south,condition_type::dirichlet,&solVal);
    bcInfo2.addCondition(0,boundary::north,condition_type::dirichlet,&solVal);
    bcInfo2.addCondition(0,boundary::west,condition_type::dirichlet,&solVal);
    bcInfo2.addCondition(0,boundary::east,condition_type::dirichlet,&solVal);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    // creating assembler
    gsElPoissonAssembler<real_t> assembler(geometry,basis,bcInfo,source);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLU solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::LU solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLU solver in " << clock.stop() <<"s.\n";
#endif

    gsMultiPatch<> solutionLaplace;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),solutionLaplace);


    // source function, rhs
    gsConstantFunction<> g(0.,2);
    gsElPoissonAssembler<real_t> assembler2(geometry,basis,bcInfo2,g);
    gsInfo<<"Assembling...\n";
    clock.restart();
    assembler2.assemble();
    gsInfo << "Assembled a system with "
           << assembler2.numDofs() << " dofs in " << clock.stop() << "s.\n";
    //assembler2.assembleRhs(solutionLaplace);

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    solVector = solver.solve(assembler2.rhs());
    gsInfo << "Solved the system with PardisoLU solver in " << clock.stop() <<"s.\n";
#else
    solVector = solver.solve(assembler2.rhs());
    gsInfo << "Solved the system with EigenLU solver in " << clock.stop() <<"s.\n";
#endif

    gsMultiPatch<> solution;
    assembler2.constructSolution(solVector,assembler2.allFixedDofs(),solution);

    //=============================================//
                  // Output //
    //=============================================//

    if (numPlotPoints > 0) // visualization
    {
        // constructing isogeometric field (geometry + solution)
        gsField<> mainField(geometry,solution);
        gsField<> auxField(geometry,solutionLaplace);
        gsField<> mainAnalytical(geometry,solVal,false);
        gsField<> auxAnalytical(geometry,laplace,false);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Main"] = &mainField;
        fields["Auxiliary"] = &auxField;
        fields["MainAnalytical"] = &mainAnalytical;
        fields["AuxiliaryAnalytical"] = &auxAnalytical;
        gsWriteParaviewMultiPhysics(fields,"biharmonicAnn",numPlotPoints);
        gsInfo << "Open \"biharmonicAnn.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
