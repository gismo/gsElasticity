/// This is an example of using the nonlinear elasticity solver on a 3D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"terrific.xml";
    index_t numUniRef = 0; // number of h-refinements
    index_t maxNumIteration = 100;
    real_t tolerance = 1e-12;
    index_t numPlotPoints = 10000;
    index_t materialLaw = material_law::saint_venant_kirchhoff;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("i","iter","Max number of iterations for Newton's method",maxNumIteration);
    cmd.addReal("t","tol","Tolerance value of Newton's method",tolerance);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - NeoHooke_ln, 2 - NeoHooke_2",materialLaw);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> f(0.,0.,0.,3);
    // surface load, neumann BC
    gsConstantFunction<> g(15e7, -10.5e7, 0,3);

    // material parameters
    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BC are imposed separately for every component (coordinate)
    for (index_t d = 0; d < 3; d++)
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
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,f);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);
    assembler.options().setInt("MaterialLaw",materialLaw);

    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsElNewton<real_t> newton(assembler);
    newton.options().setInt("MaxIters",maxNumIteration);
    newton.options().setReal("AbsTol",tolerance);
    newton.options().setInt("Verbosity",newton_verbosity::all);


    gsMultiPatch<> solutionLinear;
    index_t iterationNumber = 0;
    newton.setPostProcessingFunction([&](const gsMatrix<> & matrix) {
        if (iterationNumber == 0)
            assembler.constructSolution(matrix,solutionLinear);

        iterationNumber++;
    });

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    newton.solve();

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<> solutionNonlinear;
    assembler.constructSolution(newton.solution(),solutionNonlinear);

    //=============================================//
                  // Output //
    //=============================================//

    gsField<> nonlinearSolutionField(geometry,solutionNonlinear);
    gsField<> linearSolutionField(geometry,solutionLinear);

    gsInfo << "Plotting the output to the Paraview file \"terrific.pvd\"...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Deformation (nonlinElast)"] = &nonlinearSolutionField;
    fields["Deformation (linElast)"] = &linearSolutionField;
    gsWriteParaviewMultiPhysics(fields,"terrific",numPlotPoints);
    gsInfo << "Done. Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
