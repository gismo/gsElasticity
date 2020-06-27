/// This is a simple numerical example of modeling a muscle fiber using the nonlinear elasticity solver
/// in a mixed displacement-pressure formulation. It corresponds to Example 5.1 from the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a muscle fiber benchmark with mixed nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/muscleBeam.xml";
    real_t youngsModulusMuscle = 3.0e5; // shear modulus 1e5;
    //real_t youngsModulusTendon = 3.0e6; // shear modulus 1e6;
    real_t poissonsRatioMuscle = 0.5;
    //real_t poissonsRatioTendon = 0.5;
    index_t numUniRefDirX = 3;
    index_t numUniRef = 1;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    index_t numPlotPoints = 64000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a muscle fiber benchmark with mixed nonlinear elasticity solver.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRef);
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
    geometry.computeTopology();

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
    for (index_t i = 0; i < numUniRefDirX; ++i)
    {
        static_cast<gsTensorNurbsBasis<3,real_t> &>(basisDisplacement.basis(0)).knots(0).uniformRefine();
        static_cast<gsTensorNurbsBasis<3,real_t> &>(basisPressure.basis(0)).knots(0).uniformRefine();
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
    for (index_t d = 0; d < 3; ++d)
    {
        bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
        bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,d);
    }
    // source function, rhs
    gsConstantFunction<> g(0.,0.,5.0e4,3);

    gsFunctionExpr<> tendonMuscle("16*(1-x)^2*x^2",3);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulusMuscle);
    assembler.options().setReal("PoissonsRatio",poissonsRatioMuscle);
    assembler.options().setInt("MaterialLaw",material_law::mixed_neo_hooke_ln);
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

    if (numPlotPoints > 0) // visualization
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(geometry,displacement);
        gsField<> pressureField(geometry,pressure);
        gsField<> muscleTendonField(geometry,tendonMuscle,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["Pressure"] = &pressureField;
        fields["Muscle/tendon"] = &muscleTendonField;
        gsWriteParaviewMultiPhysics(fields,"muscleBeam",numPlotPoints);
        gsInfo << "Open \"muscleBeam.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
