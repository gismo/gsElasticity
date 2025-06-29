/// This is a simple numerical example of modeling active muscle behavior. It is based on the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
/// Muscle behavior is modeled with the incompressible nonlinear elasticity equations (pressure-displacement formulation)
/// with tendon and muscle materials for the passive part and a fiber-based active response part.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsMuscleAssembler.h>
#include <gsElasticity/src/gsMassAssembler.h>
#include <gsElasticity/src/gsElTimeIntegrator.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a simulation of active muscle behavior.\n";

    //=====================================//
                // Input //
    //=====================================//

    real_t youngsModulusMuscle = 3.0e5; // shear modulus 1e5;
    real_t youngsModulusTendon = 3.0e6; // shear modulus 1e6;
    real_t poissonsRatioMuscle = 0.5;
    real_t poissonsRatioTendon = 0.5;
    real_t density = 9e2;
    real_t gravityAcc = 0.;
    real_t maxMuscleStress = 3.0e5; // maximum stress produced at the optimal fiber stretch
    real_t optFiberStretch = 1.3;   // optimal fiber stretch
    real_t deltaW = 0.3;            // shape parameter of the active reponse function
    real_t powerNu = 4.0;           // another shape parameter of the active reponse function
    // 1 = muscle, 0 = tendon
    gsFunctionExpr<> tendonMuscleSinglePatch("16*(1-x)^2*x^2",3);
    bool rightOrLeft = true;        // true - simulate right muscle; false simulate left muscle

    // direction of muscle fibers in the parametric domain
    gsVector<> fiberDirection(3);
    fiberDirection << 1.,0.,0.;
    // space discretization
    index_t numUniRefDirX = 0;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    // output
    index_t numPlotPoints = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a simulation of active muscle behavior.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRefDirX);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation applications",numDegElev);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgridOrTaylorHood);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("m","maxStress","Maximum stress produced at the optimal fiber stretch",maxMuscleStress);
    cmd.addSwitch("left","Simulate left biceps head (default - right)",rightOrLeft);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(rightOrLeft ? gsElasticity_DATA"/bicepsRightMP.xml" :
                               gsElasticity_DATA"/bicepsLeftMP.xml", geometry);
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
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        for (index_t i = 0; i < numUniRefDirX; ++i)
        {
            static_cast<gsTensorNurbsBasis<3,real_t> &>(basisDisplacement.basis(p)).knots(0).uniformRefine();
            static_cast<gsTensorNurbsBasis<3,real_t> &>(basisPressure.basis(p)).knots(0).uniformRefine();
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
            bcInfo.addCondition(p,rightOrLeft ? boundary::left : boundary::right,condition_type::dirichlet,nullptr,d);
    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,gravityAcc*density,3);

    gsPiecewiseFunction<> tendonMuscleDistribution;
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        tendonMuscleDistribution.addPiece(tendonMuscleSinglePatch);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // creating stiffness assembler
    gsMuscleAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,gravity,tendonMuscleDistribution,fiberDirection);
    assembler.options().setReal("MuscleYoungsModulus",youngsModulusMuscle);
    assembler.options().setReal("MusclePoissonsRatio",poissonsRatioMuscle);
    assembler.options().setReal("TendonYoungsModulus",youngsModulusTendon);
    assembler.options().setReal("TendonPoissonsRatio",poissonsRatioTendon);
    assembler.options().setReal("MaxMuscleStress",maxMuscleStress);
    assembler.options().setReal("OptFiberStretch",optFiberStretch);
    assembler.options().setReal("DeltaW",deltaW);
    assembler.options().setReal("PowerNu",powerNu);
    assembler.options().setReal("Alpha",1.);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    // setting Newton's method
    gsIterative<real_t> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setReal("RelTol",1e-4); // Should be 1e-9
    solver.options().setInt("Solver",linear_solver::LDLT);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    solver.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                    // Output //
    //=============================================//

    // displacement as an isogeometric displacement field
    gsMultiPatch<> displacement,pressure;
    assembler.constructSolution(solver.solution(),solver.allFixedDofs(),displacement,pressure);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

    if (numPlotPoints > 0) // visualization
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),displacement);
        gsField<> pressureField(assembler.patches(),pressure);
        gsField<> stressField(assembler.patches(),stresses,true);
        gsField<> muscleTendonField(geometry,tendonMuscleDistribution,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["Pressure"] = &pressureField;
        fields["von Mises"] = &stressField;
        fields["Muscle/tendon"] = &muscleTendonField;
        gsWriteParaviewMultiPhysics(fields,rightOrLeft ? "bicepsRightMS" : "bicepsLeftMS",numPlotPoints);
        gsInfo << "Open \" " << (rightOrLeft ? "bicepsRightMS" : "bicepsLeftMS")
               << ".pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
