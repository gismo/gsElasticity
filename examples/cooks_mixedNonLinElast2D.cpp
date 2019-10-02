/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver with a mixed displacement-pressure formulation.
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

    gsInfo << "This is Cook's membrane benchmark with mixed nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cooks.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 1; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t youngsModulus = 240.565e6;
    real_t poissonsRatio = 0.4;
    bool subgrid = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("s","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

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
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);

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

    // displacement as an isogeometric displacement field
    gsMultiPatch<> displacement,pressure;
    assembler.constructSolution(newton.solution(),displacement,pressure);

    //=============================================//
                  // Visualization //
    //=============================================//

    if (numPlotPoints)
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),displacement);
        gsField<> pressureField(assembler.patches(),pressure);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["Pressure"] = &pressureField;
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
