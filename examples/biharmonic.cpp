/// This is an example of using the biharmonic equation solver.
///
/// Authors: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsBiharmonicAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the biharmonic solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"unitSquare.xml";
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
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basisMain(geometry);
    gsMultiBasis<> basisAux(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisMain.degreeElevate();
        basisAux.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisMain.uniformRefine();
        basisAux.uniformRefine();
    }
    // additional velocity refinement for stable mixed FEM
    //if (!subgridOrTaylorHood) // subgrid
    //    basisMain.uniformRefine();
    //else // Taylor-Hood
    //    basisAux.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> f(-1.,2);
    // surface load, neumann BC

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0);
    //bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,0,1);
    //bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
    //bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    //bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    // creating assembler
    gsBiharmonicAssembler<real_t> assembler(geometry,basisMain,basisAux,bcInfo,f);
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

    //=============================================//
                  // Output //
    //=============================================//

    // constructing solution as an IGA function
    gsMultiPatch<> solutionMain, solutionAux;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),solutionMain,solutionAux);

    if (numPlotPoints > 0) // visualization
    {
        // constructing isogeometric field (geometry + solution)
        gsField<> mainField(assembler.patches(),solutionMain);
        gsField<> auxField(assembler.patches(),solutionAux);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Main"] = &mainField;
        fields["Auxiliary"] = &auxField;
        gsWriteParaviewMultiPhysics(fields,"biharmonic",numPlotPoints);
        gsInfo << "Open \"biharmonic.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
