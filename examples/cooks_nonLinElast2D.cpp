/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cooks.xml";
    index_t numUniRef = 4; // number of h-refinements
    index_t numKRef = 1; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t youngsModulus = 240.565e6;
    real_t poissonsRatio = 0.4;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("s","point","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisDisplacement.degreeElevate();
        basisDisplacement.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
    }

    // neumann BC
    gsConstantFunction<> f(0.,625e4,2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (index_t d = 0; d < 2; ++d)
        bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
    bcInfo.addCondition(0,boundary::east,condition_type::neumann,&f);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsNewton<real_t> newton(assembler);
    newton.options().setInt("Verbosity",newton_verbosity::all);
    newton.options().setInt("Solver",linear_solver::LDLT);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<> displacement;
    assembler.constructSolution(newton.solution(),displacement);

    //=============================================//
                  // Output //
    //=============================================//

    if (numPlotPoints)
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),displacement);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        gsWriteParaviewMultiPhysics(fields,"cooks",numPlotPoints);
        gsInfo << "Open \"cooks.pvd\" in Paraview for visualization.\n";
    }

    //=============================================//
                  // Validation //
    //=============================================//

    gsMatrix<> A(2,1);
    A << 1.,1.;
    gsInfo << "Y-displacement of the top-right corner: " << displacement.patch(0).eval(A)(1,0) << std::endl;

    return 0;
}
