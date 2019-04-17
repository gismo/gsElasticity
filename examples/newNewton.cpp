/// This is an example of using the nonlinear elasticity solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityNewton2.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/lshape.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numDegElevate = 1; // number of p-refinements
    index_t maxNumIteration = 100;
    real_t tolerance = 1e-12;
    index_t numPlotPoints = 10000;
    index_t materialLaw = material_law::saint_venant_kirchhoff;
    index_t numIncSteps = 1;
    index_t save = newtonSave::onlyFinal;
    index_t verbosity = newtonVerbosity::all;
    bool plotDeform = true;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the nonlinear elasticity solver in 2D.");
    cmd.addInt("i","iter","Max number of iterations for Newton's method",maxNumIteration);
    cmd.addReal("t","tol","Tolerance value of Newton's method",tolerance);
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - NeoHooke_ln, 2 - NeoHooke_2",materialLaw);
    cmd.addInt("n","num","Number incremental step",numIncSteps);
    cmd.addInt("s","save","Save",save);
    cmd.addInt("v","verbosity","Verbosity",verbosity);
    cmd.addSwitch("p","plot","Plot",plotDeform);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // boundary displacement in y-direction, dirichlet BC
    gsConstantFunction<> f(0.,1e9,2);

    // material parameters
    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0); // Dirichlet BC are defined componentwise
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1); // Last number is a component (coordinate) index
    bcInfo.addCondition(2,boundary::north,condition_type::neumann,&f);

    //=============================================//
                  // Setting //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basis(geometry);
    for (index_t i = 0; i < numDegElevate; ++i)
        basis.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);
    assembler.options().setInt("MaterialLaw",materialLaw);

    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsElasticityNewton2<real_t> newton(assembler);
    newton.options().setInt("NumIncStep",numIncSteps);
    newton.options().setInt("Save",save);
    newton.options().setInt("Verbosity",verbosity);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    newton.solve();
    if (plotDeform)
    {
        for (index_t i = 0; i < numDegElevate; ++i)
            geometry.degreeElevate();
        for (index_t i = 0; i < numUniRef; ++i)
            geometry.uniformRefine();
        newton.plotDeformation(geometry,"lshapeDeform",0);
    }

    // solution to the nonlinear problem as an isogeometric displacement field
    const gsMultiPatch<> solutionNonlinear = newton.solution();
    // solution to the linear problem as an isogeometric displacement field
    //const gsMultiPatch<> solutionLinear = newton.linearSolution();

    //=============================================//
                  // Output //
    //=============================================//

    gsField<> nonlinearSolutionField(geometry,solutionNonlinear);
    //gsField<> linearSolutionField(geometry,solutionLinear);

    gsInfo << "Plotting the output to the Paraview file \"lshape.pvd\"...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Deformation (nonlinElast"] = &nonlinearSolutionField;
    //fields["Deformation (linElast)"] = &linearSolutionField;
    gsWriteParaviewMultiPhysics(fields,"lshape",numPlotPoints);
    gsInfo << "Done. Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
