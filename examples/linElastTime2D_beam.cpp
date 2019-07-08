/// This is an example of using the time-dependent linear elasticity solver on a 2D geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElMassAssembler.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the time-dependent linear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/unitSquare.xml";
    index_t numUniRef = 0; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.1,2);

    // material parameters
    real_t youngsModulus = 200.;//74e9;
    real_t poissonsRatio = 0.33;
    real_t density = 1.;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);

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

    // creating assemblers
    gsElasticityAssembler<real_t> stiffAssembler(geometry,basis,bcInfo,g);
    stiffAssembler.options().setReal("YoungsModulus",youngsModulus);
    stiffAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    stiffAssembler.options().setInt("DirichletValues",dirichlet::interpolation);

    gsElMassAssembler<real_t> massAssembler(geometry,basis,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    gsInfo<<"Assembling stiffness...\n";
    gsStopwatch clock;
    clock.restart();
    stiffAssembler.assemble();
    gsInfo << "Assembled the stiffness matrix and the rhs vector with "
           << stiffAssembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo<<"Assembling mass...\n";
    clock.restart();
    massAssembler.assemble();
    gsInfo << "Assembled the mass matrix with "
           << massAssembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << massAssembler.matrix().toDense() << std::endl;
    return 0;
}
