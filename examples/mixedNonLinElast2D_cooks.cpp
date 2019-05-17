/// This is an example of using the mixed nonlinear elasticity solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityNewtonDeLuxe.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the mixed nonlinear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cooks.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 1; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t poissonsRatio = 0.4;
    index_t numSteps = 1;
    bool subgrid = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addInt("i","iter","Number of incremental loading steps",numSteps);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // neumann BC
    gsConstantFunction<> f(0.,625e4,2);

    // material parameters
    real_t youngsModulus = 240.565e6;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    bcInfo.addCondition(0,boundary::east,condition_type::neumann,&f);

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisDisplacement.degreeElevate();
        basisPressure.degreeElevate();
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    if (subgrid)
        basisDisplacement.uniformRefine();
    else
        basisDisplacement.degreeElevate();

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    assembler.options().setInt("MaterialLaw",1);

    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsElasticityNewtonDeLuxe<real_t> newton(assembler);
    newton.options().setInt("Verbosity",newtonVerbosity2::all);
    newton.options().setInt("Save",newtonSave2::firstAndLastPerIncStep);
    newton.options().setInt("NumIncStep",numSteps);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    newton.solve();

    // displacement as an isogeometric displacement field
    const gsMultiPatch<> & displacement = newton.displacement();
    // pressure as an isogeometric displacement field
    const gsMultiPatch<> & pressure = newton.pressure();

    //=============================================//
                  // Output //
    //=============================================//

    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(),displacement);
    gsField<> pressureField(assembler.patches(),pressure);

    gsInfo << "Plotting the output to the Paraview file \"cooks.pvd\"...\n";
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &displacementField;
    fields["Pressure"] = &pressureField;
    gsWriteParaviewMultiPhysics(fields,"cooks",numPlotPoints);
    gsInfo << "Done. Use Warp-by-Vector filter in Paraview to deform the geometry.\n";

    return 0;
}
