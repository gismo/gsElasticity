/// This is a simple numerical example of modeling a muscle fiber using the nonlinear elasticity solver
/// in a mixed displacement-pressure formulation. It corresponds to Example 5.1 from the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
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
    real_t density = 9e2;
    real_t loading = 5.0e4;
    // space discretization
    index_t numUniRefDirX = 3;
    index_t numUniRef = 1;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    // time integration
    real_t timeSpan = 1;
    real_t timeStep = 0.1;
    // output
    index_t numPlotPoints = 64000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a muscle fiber benchmark with mixed nonlinear elasticity solver.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRef);
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
    gsConstantFunction<> gravity(0.,0.,loading,3);

    gsFunctionExpr<> tendonMuscle("16*(1-x)^2*x^2",3);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    gsInfo << "Solvers...\n";
    // creating stiffness assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulusMuscle);
    assembler.options().setReal("PoissonsRatio",poissonsRatioMuscle);
    assembler.options().setInt("MaterialLaw",material_law::mixed_neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisDisplacement,bcInfo,gravity);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsElTimeIntegrator<real_t> timeSolver(assembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    timeSolver.options().setInt("Verbosity",solver_verbosity::all);

    //=============================================//
            // Setting output & auxilary//
    //=============================================//

    gsInfo << "Output...\n";
    // displacement field
    gsMultiPatch<> displacement;
    // stress field
    gsPiecewiseFunction<> stresses;
    // constructing an IGA field (geometry + solution)
    gsField<> dispField(geometry,displacement);
    //gsField<> stressField(assembler.patches(),stresses,true);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;
    //fields["von Mises"] = &stressField;

    //std::ofstream logFile;
    //logFile.open("muscleBeam.txt");
    //logFile << "# simTime dispAx dispAy compTime numIters\n";

    gsProgressBar bar;
    gsStopwatch totalClock;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    gsInfo << "Init conod...\n";
    // set initial conditions
    gsInfo << "Set disp...\n";
    timeSolver.setDisplacementVector(gsMatrix<>::Zero(assembler.numDofs(),1));
    gsInfo << "Set vel...\n";
    timeSolver.setVelocityVector(gsMatrix<>::Zero(assembler.numDofs(),1));

    gsInfo << "Constr...\n";
    assembler.constructSolution(timeSolver.displacementVector(),timeSolver.allFixedDofs(),displacement);
    //assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);
    //writeLog(logFile,displacement,0.,0.,0);
    // plotting initial displacement

    gsInfo << "Plot...\n";
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
        assembler.constructSolution(timeSolver.displacementVector(),timeSolver.allFixedDofs(),displacement);

        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"muscleBeam",collection,i+1,numPlotPoints);
        //writeLog(logFile,assembler,velocity,pressure,simTime,compTime,timeSolver.numberIterations());
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

    //logFile.close();
    //gsInfo << "Log file created in \"flappingBeam_CSM3.txt\".\n";


    return 0;
}
