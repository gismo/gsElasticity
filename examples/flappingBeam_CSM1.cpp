/// This is the structural solver benchmark CSM1 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Benchmark CSM1: steady-state deformation of an elastic beam.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 1; // number of k-refinements
    index_t numPlotPoints = 10000;
    real_t poissonsRatio = 0.4;
    real_t youngsModulus = 1.4e6;
    real_t density = 1.0e3;
    real_t gravitationalAcc = 2.0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CSM1: steady-state deformation of an elastic beam.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
                  // Setting solver //
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
        basisDisplacement.uniformRefine();

    // boundary conditions
    gsBoundaryConditions<> bcInfo; // numbers are: patch, function pointer (nullptr) for displacement, displacement component
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);

    // gravity, rhs
    gsConstantFunction<> gravity(0.,-1*gravitationalAcc*density,2);

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    // setting Newton's method
    gsNewton<real_t> newton(assembler);
    newton.options().setInt("Verbosity",newton_verbosity::all);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    // solution as an isogeometric displacement field
    gsMultiPatch<> solution;
    assembler.constructSolution(newton.solution(),solution);

    //=============================================//
                  // Validation //
    //=============================================//

    gsMatrix<> A(2,1);
    A << 1,0.5;
    gsInfo << "Displacement of the point A:\n" << solution.patch(0).eval(A) << std::endl;

    //=============================================//
                  // Visualization //
    //=============================================//

    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(),solution);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &displacementField;
    gsWriteParaviewMultiPhysics(fields,"flappingBeam_CSM1",numPlotPoints);
    gsInfo << "Open \"flappingBeam_CSM1.pvd\" in Paraview for visualization.\n";

    return 0;
}
