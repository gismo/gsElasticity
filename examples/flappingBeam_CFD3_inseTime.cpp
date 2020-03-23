/// This is the incompressible Navier-Stokes solver benchmark CFD3 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
///
/// INSE can be solved the Newton or the semi-implicit (a.k.a. IMEX) time integration schemes
/// (see the book by Volker John: Finite Element Methods for Incompressible Flow Problems, Springer, 2016).
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
///
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

void refineBoundaryLayer(gsMultiBasis<> & velocity, gsMultiBasis<> & pressure)
{
    gsMatrix<> boxEast(2,2);
    boxEast << 0.8,1.,0.,0.;
    velocity.refine(0,boxEast);
    pressure.refine(0,boxEast);

    gsMatrix<> boxSouth(2,2);
    boxSouth << 0.,0.,0.,0.2;
    velocity.refine(1,boxSouth);
    velocity.refine(3,boxSouth);
    pressure.refine(1,boxSouth);
    pressure.refine(3,boxSouth);

    gsMatrix<> boxNorth(2,2);
    boxNorth << 0.,0.,0.8,1.;
    velocity.refine(2,boxNorth);
    velocity.refine(4,boxNorth);
    pressure.refine(2,boxNorth);
    pressure.refine(4,boxNorth);

    gsMatrix<> boxWest(2,2);
    boxWest << 0.,0.2,0.,0.;
    velocity.refine(5,boxWest);
    pressure.refine(5,boxWest);
}

void writeLog(std::ofstream & ofs, const gsNsAssembler<real_t> & assembler,
              const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure,
              real_t simTime, real_t compTime, index_t numIters)
{
    // computing force acting on the surface of the submerged structure
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);

    // compute the pressure difference between the front and the end points of the structure
    gsMatrix<> A(2,1);
    A << 1.,0.5;
    A = pressure.patch(0).eval(A);
    gsMatrix<> B(2,1);
    B << 0.,0.5;
    B = pressure.patch(5).eval(B);

    // print: simTime drag lift pressureDiff compTime numIters
    ofs << simTime << " " << force.at(0) << " " << force.at(1) << " "
        << A.at(0)-B.at(0)<< " " << compTime << " " << numIters << std::endl;
}

int main(int argc, char* argv[]){
    gsInfo << "Benchmark CFD3: transient flow of an incompressible fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flappingBeam_flow.xml";
    real_t viscosity = 0.001; // kinematic viscosity
    real_t meanVelocity = 2; // inflow velocity
    real_t density = 1.0e3;
    // space discretization
    index_t numUniRef = 3;
    index_t numDegElev = 0;
    index_t numBLRef = 1;
    bool subgridOrTaylorHood = false;
    // time integration
    real_t timeSpan = 10;
    real_t timeStep = 0.01;
    real_t theta = 0.5;
    bool imexOrNewton = false;
    bool warmUp = false;
    // output
    index_t numPlotPoints = 900;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CFD3: transient flow of an incompressible fluid.");
    cmd.addInt("r","refine","Number of uniform h-refinements",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevations",numDegElev);
    cmd.addInt("l","blayer","Number of additional h-refinements in the boundary layer",numBLRef);
    cmd.addSwitch("e","element","Mixed element: false = subgrid (default), true = Taylor-Hood",subgridOrTaylorHood);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addReal("f","theta","Time integration parameter: 0 - exp.Euler, 1 - imp.Euler, 0.5 - Crank-Nicolson",theta);
    cmd.addSwitch("i","intergration","Time integration scheme: false = IMEX (default), true = Newton",imexOrNewton);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of sampling points per patch for Paraview (0 = no plotting)",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements with the "
           << (imexOrNewton ? "Newton " : "IMEX ") << "time integration scheme.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisVelocity(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    if (!subgridOrTaylorHood)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0); // x-component
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);       // y-component
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(5,boundary::west,condition_type::dirichlet,0,d);
        bcInfo.addCondition(6,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(6,boundary::north,condition_type::dirichlet,0,d);
    }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // creating stiffness assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setReal("Density",density);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsNsTimeIntegrator<real_t> timeSolver(assembler,massAssembler);
    timeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    timeSolver.options().setReal("Theta",theta);

    //=============================================//
            // Setting output & auxilary//
    //=============================================//

    // solution as two isogeometric fields
    gsMultiPatch<> velocity, pressure;
    // isogeometric fields (geometry + solution)
    gsField<> velocityField(geometry,velocity);
    gsField<> pressureField(geometry,pressure);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;
    // paraview collection of time steps
    gsParaviewCollection collection("flappingBeam_CFD3");

    std::ofstream logFile;
    logFile.open("flappingBeam_CFD3.txt");

    gsProgressBar bar;
    gsStopwatch iterClock, totalClock;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    // I change the Dirichlet DoFs for warming up, so I save them here for later
    gsMatrix<> inflowDDoFs;
    assembler.getFixedDofs(0,boundary::west,inflowDDoFs);
    assembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    timeSolver.setSolutionVector(gsMatrix<>::Zero(assembler.numDofs(),1));
    timeSolver.setFixedDofs(assembler.allFixedDofs());
    timeSolver.initialize();

    // consruct and plot initial velocity
    assembler.constructSolution(timeSolver.solutionVector(),timeSolver.allFixedDofs(),velocity,pressure);
    writeLog(logFile,assembler,velocity,pressure,0.,0.,0);
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CFD3",collection,0,numPlotPoints);

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

        // change time step for the initial warm-up phase
        real_t tStep = (warmUp && simTime < 2.) ? 0.1 : timeStep;
        // smoothly change the inflow boundary condition
        if (simTime < 2.)
            assembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*(simTime+tStep)/2.))/2);

        timeSolver.makeTimeStep(tStep);
        // construct solution; timeSolver already knows the new Dirichlet BC
        assembler.constructSolution(timeSolver.solutionVector(),timeSolver.allFixedDofs(),velocity,pressure);

        simTime += tStep;
        numTimeStep++;
        compTime += iterClock.stop();

        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CFD3",collection,numTimeStep,numPlotPoints);
        writeLog(logFile,assembler,velocity,pressure,simTime,compTime,timeSolver.numberIterations());
    }

    //=============================================//
                // Final touches //
    //=============================================//

    gsInfo << "Simulation time: " + secToHMS(compTime) << " (total time: " + secToHMS(totalClock.stop()) + ")\n";

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"flappingBeam_CFD3.pvd\" in Paraview for visualization.\n";
    }
    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_CFD3.txt\".\n";

    return 0;
}
