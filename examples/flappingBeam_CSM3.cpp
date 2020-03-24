/// This is the structural solver benchmark CSM3 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

void writeLog(std::ofstream & file, const gsMultiPatch<> & displacement,
              real_t simTime, real_t compTime, real_t numIters)
{
    // evaluating displacement at the point A
    gsMatrix<> point(2,1);
    point << 1., 0.5;
    gsMatrix<> dispA = displacement.patch(0).eval(point);

    // print: simTime dispAx dispAy compTime numIters
    file << simTime << " " << dispA.at(0) << " " << dispA.at(1) << " "
         << compTime << " " << numIters << std::endl;
}

int main(int argc, char* argv[]){
    gsInfo << "Benchmark CSM3: dynamic deformation of an elastic beam.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    real_t poissonsRatio = 0.4;
    real_t youngsModulus = 1.4e6;
    real_t density = 1.0e3;
    real_t gravitationalAcc = 2.0;
    // space discretization
    index_t numUniRef = 3;
    index_t numDegElev = 0;
    // time integration
    real_t timeSpan = 2;
    real_t timeStep = 0.01;
    // output
    index_t numPlotPoints = 1000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CSM3: dynamic deformation of an elastic beam.");
    cmd.addReal("g","graviry","Gravitational acceleration",gravitationalAcc);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addInt("p","points","Number of sampling points to plot to Paraview",numPlotPoints);
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
    gsBoundaryConditions<> bcInfo; // numbers are: patch, function pointer (nullptr) for displacement, displacement component
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);

    // gravity, rhs
    gsConstantFunction<> gravity(0.,gravitationalAcc*density,2);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // creating stiffness assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,bcInfo,gravity);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
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
    gsMultiPatch<> displacement;

    // constructing an IGA field (geometry + solution)
    gsField<> dispField(geometry,displacement);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;

    std::ofstream logFile;
    logFile.open("flappingBeam_CSM3.txt");
    logFile << "# simTime dispAx dispAy compTime numIters\n";

    gsProgressBar bar;
    gsStopwatch iterClock, totalClock;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    // set initial conditions
    timeSolver.setDisplacementVector(gsMatrix<>::Zero(assembler.numDofs(),1));
    timeSolver.setVelocityVector(gsMatrix<>::Zero(assembler.numDofs(),1));

    assembler.constructSolution(timeSolver.displacementVector(),displacement);

    // plotting initial displacement
    gsParaviewCollection collection("flappingBeam_CSM3");
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CSM3",collection,0,numPlotPoints);


    //=============================================//
                  // Solving //
    //=============================================//

    real_t simTime = 0.;
    real_t numTimeStep = 0;
    real_t compTime = 0.;

    gsInfo << "Running the simulation...\n";
    totalClock.restart();
    while (simTime < timeSpan)
    {
        bar.display(simTime/timeSpan);
        iterClock.restart();

        timeSolver.makeTimeStep(timeStep);
        assembler.constructSolution(timeSolver.displacementVector(),displacement);

        compTime += iterClock.stop();
        simTime += timeStep;
        numTimeStep++;

        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CSM3",collection,numTimeStep,numPlotPoints);
        writeLog(logFile,displacement,simTime,compTime,timeSolver.numberIterations());
    }

    //=============================================//
                // Final touches //
    //=============================================//

    gsInfo << "Simulation time: " + secToHMS(compTime) << " (total time: " + secToHMS(totalClock.stop()) + ")\n";

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"flappingBeam_CSM3.pvd\" in Paraview for visualization.\n";
    }

    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_CSM3.txt\".\n";

    return 0;
}
