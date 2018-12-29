/// This is an example of using the nonlinear elasticity solver on a 3D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityNewton.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"terrific.xml";
    int numUniRef = 0; // number of h-refinements
    int numDegElevate = 0; // number of p-refinements
    int maxNumIteration = 100;
    real_t tolerance = 1e-12;
    int numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","prefine","Number of degree elevation application",numDegElevate);
    cmd.addInt("i","iter","Max number of iterations for Newton's method",maxNumIteration);
    cmd.addReal("t","tol","Tolerance value of Newton's method",tolerance);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> f(0.,0.,0.,3);
    // surface load, neumann BC
    gsConstantFunction<> g(20e7, -14e7, 0,3);

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

    // setting Newton's method
    gsElasticityNewton<real_t> newton(assembler);
    newton.setMaxIterations(maxNumIteration);
    newton.setTolerance(tolerance);
    newton.setVerbosity(true);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    newton.solve();

    // solution to the nonlinear problem as an isogeometric displacement field
    const gsMultiPatch<> solutionNonlinear = newton.solution();
    // solution to the linear problem as an isogeometric displacement field
    const gsMultiPatch<> solutionLinear = newton.linearSolution();

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
