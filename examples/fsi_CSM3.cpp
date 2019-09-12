/// This is the structural solver benchmark CSM3 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

void validation(std::ofstream & file, const gsMultiPatch<> & displacement)
{
    // evaluating displacement at the point A
    gsMatrix<> point(2,1);
    point << 1., 0.5;
    gsMatrix<> dispA = displacement.patch(0).eval(point);
    file << dispA.at(0) << " " << dispA.at(1) << std::endl;
}

int main(int argc, char* argv[]){

    gsInfo << "Benchmark CSM3: dynamic deformation of an elastic beam.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/fsi_flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 1; // number of k-refinements
    real_t poissonsRatio = 0.4;
    real_t youngsModulus = 1.4e6;
    real_t density = 1.0e3;
    real_t gravitationalAcc = 2.0;
    // roughly, period = 2 * pi * sqrt(density)/width(0.02)/sqrt(YoungsModulus) = 2.22
    real_t timeSpan = 2;
    real_t timeStep = 0.01;
    index_t numPlotPoints = 10000;
    bool validate = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CSM1: steady-state deformation of an elastic beam.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("p","points","Number of sampling points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addSwitch("x","validate","Save displacement of the point A to a text file for further analysis",validate);
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
    // set initial conditions
    timeSolver.setInitialDisplacement(gsMatrix<>::Zero(assembler.numDofs(),1));
    timeSolver.setInitialVelocity(gsMatrix<>::Zero(assembler.numDofs(),1));
    timeSolver.initialize();

    //=============================================//
            // Setting output & auxilary//
    //=============================================//

    // initial displacement field
    gsMultiPatch<> displacement;
    assembler.constructSolution(gsMatrix<>::Zero(assembler.numDofs(),1),displacement);

    // constructing an IGA field (geometry + solution)
    gsField<> dispField(geometry,displacement);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;

    // plotting initial displacement
    gsParaviewCollection collection("fsi_CSM3");
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"fsi_CSM3",collection,0,numPlotPoints);

    gsProgressBar bar;
    gsStopwatch clock;

    std::ofstream file;
    if (validate)
        file.open("fsi_CSM3.txt");

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Running the transient simulation...\n";
    clock.restart();
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));
        timeSolver.makeTimeStep(timeStep);
        assembler.constructSolution(timeSolver.displacementVector(),displacement);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"fsi_CSM3",collection,i+1,numPlotPoints);
        if (validate)
            validation(file,displacement);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"fsi_CSM3.pvd\" in Paraview for visualization.\n";
    }
    if (validate)
    {
        file.close();
        gsInfo << "Displacement of the point A over time is saved to \"fsi_CSM3.txt\".\n";
    }

    return 0;
}
