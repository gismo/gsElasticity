/// This is a simple numerical example of modeling a dynamic muscle fiber using the nonlinear elasticity solver
/// in a mixed displacement-pressure formulation. It rougly corresponds to Example 5.1 from the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsMassAssembler.h>
#include <gsElasticity/src/gsElTimeIntegrator.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a muscle fiber benchmark with a time-dependent mixed nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = gsElasticity_DATA"/muscleBeamMP.xml";
    real_t youngsModulus = 3.0e5; // shear modulus 1e5;
    real_t poissonsRatio = 0.5;
    real_t density = 9e2;
    real_t gravityAcc = -9.8;
    // space discretization
    index_t numUniRefDirX = 3;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    // time integration
    real_t timeSpan = 0.01; // was 1.0
    real_t timeStep = 0.01;
    // output
    index_t numPlotPoints = 1000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a muscle fiber benchmark with mixed nonlinear elasticity solver.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRefDirX);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation applications",numDegElev);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgridOrTaylorHood);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
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
          // Setting assemblers and solvers //
    //=============================================//

    // creating stiffness assembler for the displacement-pressure formulation
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",material_law::mixed_neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisDisplacement,bcInfo,gravity);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsElTimeIntegrator<real_t> timeSolver(assembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    timeSolver.options().setInt("Verbosity",solver_verbosity::none);

    //=============================================//
            // Setting output & auxilary//
    //=============================================//

    // displacement field
    gsMultiPatch<> displacement, pressure;
    // stress field
    gsPiecewiseFunction<> stresses;
    // constructing an IGA field (geometry + solution)
    gsField<> dispField(geometry,displacement);
    gsField<> presField(geometry,pressure);
    gsField<> stressField(assembler.patches(),stresses,true);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;
    fields["Pressure"] = &presField;
    fields["von Mises"] = &stressField;

    gsProgressBar bar;
    gsStopwatch totalClock;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    // set initial conditions
    timeSolver.setDisplacementVector(gsMatrix<>::Zero(massAssembler.numDofs(),1));
    timeSolver.setVelocityVector(gsMatrix<>::Zero(massAssembler.numDofs(),1));

    timeSolver.constructSolution(displacement,pressure);
    assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

    // plotting initial displacement
    gsParaviewCollection collection("muscleBeam");
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"muscleBeam",collection,0,numPlotPoints);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Running the simulation...\n";
    totalClock.restart();
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));

        timeSolver.makeTimeStep(timeStep);
        // construct solution; timeSolver already knows the new Dirichlet BC
        timeSolver.constructSolution(displacement,pressure);
        assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"muscleBeam",collection,i+1,numPlotPoints);
    }

    //=============================================//
                // Final touches //
    //=============================================//

    gsInfo << "Simulation time: " + secToHMS(totalClock.stop()) + "\n";

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"muscleBeam.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
