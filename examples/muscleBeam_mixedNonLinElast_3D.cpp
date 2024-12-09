/// This is a simple numerical example of modeling a muscle fiber using the nonlinear elasticity solver
/// in a mixed displacement-pressure formulation. It roughly corresponds to Example 5.1 from the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a muscle fiber benchmark with a mixed nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/muscleBeamMP.xml";
    real_t youngsModulus = 3.0e5; // shear modulus 1e5;
    real_t poissonsRatio = 0.5;
    real_t density = 9e2;
    real_t gravityAcc = -9.8;
    // spatial discretization
    index_t numUniRefDirX = 4;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    // output
    index_t numPlotPoints = 1000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a muscle fiber benchmark with mixed nonlinear elasticity solver.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRefDirX);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation applications",numDegElev);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgridOrTaylorHood);
    cmd.addInt("s","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisDisplacement.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        for (index_t i = 0; i < numUniRefDirX; ++i)
        {
            basisDisplacement.basis(p).uniformRefine(1,1,0);
            basisPressure.basis(p).uniformRefine(1,1,0);
        }
    // additional displacement refinement for stable mixed FEM
    if (!subgridOrTaylorHood) // subgrid
        basisDisplacement.uniformRefine();
    else  // Taylor-Hood
        basisDisplacement.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        for (index_t d = 0; d < 3; ++d)
        {
            bcInfo.addCondition(p,boundary::west,condition_type::dirichlet,nullptr,d);
            bcInfo.addCondition(p,boundary::east,condition_type::dirichlet,nullptr,d);
        }
    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,gravityAcc*density,3);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler for the displacement-pressure formulation
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
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

    // displacement and pressure as isogeometric fields
    gsMultiPatch<> displacement,pressure;
    assembler.constructSolution(solver.solution(),solver.allFixedDofs(),displacement,pressure);
    // construct stress field
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

    if (numPlotPoints > 0) // visualization
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(geometry,displacement);
        gsField<> pressureField(geometry,pressure);
        gsField<> stressField(assembler.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["Pressure"] = &pressureField;
        fields["von Mises"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"muscleBeam",numPlotPoints);
        gsInfo << "Open \"muscleBeam.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
