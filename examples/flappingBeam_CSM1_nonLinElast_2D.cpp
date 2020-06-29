/// This is the structural solver benchmark CSM1 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Benchmark CSM1: stationary deflection of an elastic beam.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    real_t poissonsRatio = 0.4;
    real_t youngsModulus = 1.4e6;
    real_t density = 1.0e3;
    index_t materialLaw = material_law::saint_venant_kirchhoff;
    real_t loading = 2.;
    index_t numUniRef = 3;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CSM1: stationary deflection of an elastic beam.");
    cmd.addReal("l","load","Gravitational loading acting on the beam",loading);
    cmd.addInt("m","matlaw","Material law: 0 - St.V.-K., 1 - neoHookeLn, 2 - neoHookeQuad",materialLaw);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);

    // gravity, rhs
    gsConstantFunction<> gravity(0.,loading*density,2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",materialLaw);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative<real_t> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setInt("Solver",linear_solver::LDLT);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    solver.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";  

    //=============================================//
                      // Output //
    //=============================================//

    // solution as an isogeometric displacement field
    gsMultiPatch<> solution;
    assembler.constructSolution(solver.solution(),solver.allFixedDofs(),solution);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_components::von_mises);

    if (numPlotPoints > 0) // visualization
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),solution);
        gsField<> stressField(assembler.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["von Mises"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"flappingBeam_CSM1",numPlotPoints);
        gsInfo << "Open \"flappingBeam_CSM1.pvd\" in Paraview for visualization.\n";
    }

    // validation
    gsMatrix<> A(2,1);
    A << 1.,0.5;
    A = solution.patch(0).eval(A);
    gsInfo << "X-displacement of the point A: " << A.at(0) << std::endl;
    gsInfo << "Y-displacement of the point A: " << A.at(1) << std::endl;

    return 0;
}
