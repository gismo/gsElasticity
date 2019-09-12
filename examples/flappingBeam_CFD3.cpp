/// This is the incompressible Navier-Stokes solver benchmark CFD3 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006
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

void validation(std::ofstream & ofs, const gsNsAssembler<real_t> & assembler, real_t time,
                const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure)
{
    // computing force acting on the surface of the cylinder
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);
    // print time-drag-lift
    ofs << time << " " << force.at(0) << " " << force.at(1) << std::endl;
}

int main(int argc, char* argv[]){
    gsInfo << "Benchmark CFD3: transient flow of an incompressible fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    real_t viscosity = 0.001;
    real_t meanVelocity = 2;
    real_t density = 1.0e3;
    bool subgrid = false;
    bool supg = false;
    real_t timeSpan = 10;
    real_t timeStep = 0.01;
    index_t numPlotPoints = 900;
    bool validate = true;
    real_t warmUpTimeSpan = 2.;
    real_t warmUpTimeStep = 0.1;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Benchmark CFD3: transient flow of an incompressible fluid.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Use SUPG stabilization (testing)",supg);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addSwitch("x","validate","Save lift and drag over time to a text file for further analysis",validate);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
                  // Setting solver //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisVelocity(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    // additional velocity refinement for stable mixed FEM
    if (subgrid)
        basisVelocity.uniformRefine(); // subgrid
    else
        basisVelocity.degreeElevate(); // Taylor-Hood

    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0); // x-component
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1); // y-component
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

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setReal("Density",density);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    assembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsNsTimeIntegrator<real_t> timeSolver(assembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration_NS::implicit_euler_newton);


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


    std::ofstream file;
    if (validate)
        file.open("flappingBeam_CFD3.txt");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
                   // Warming up //
    //=============================================//

    assembler.constructSolution(gsMatrix<>::Zero(assembler.numDofs(),1),velocity,pressure);
    gsMatrix<> inflowDDoFs = velocity.patch(0).boundary(boundary::west)->coefs();

    assembler.setDirichletDofs(0,boundary::west,inflowDDoFs*0.);
    assembler.constructSolution(gsMatrix<>::Zero(assembler.numDofs(),1),velocity,pressure);
    // plotting initial condition
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CFD3",collection,0,numPlotPoints);

    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    gsMatrix<> solVector = gsMatrix<>::Zero(assembler.numDofs(),1);
    for (index_t i = 0; i < index_t(warmUpTimeSpan/warmUpTimeStep); ++i)
    {
        bar.display(i+1,index_t(warmUpTimeSpan/warmUpTimeStep));
        assembler.setDirichletDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*warmUpTimeStep*(i+1)/warmUpTimeSpan))/2);
        solVector = timeSolver.oseenFSI(velocity,warmUpTimeStep,solVector);
        assembler.constructSolution(solVector,velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CFD3",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    // set initial conditions
    timeSolver.setInitialSolution(solVector);
    timeSolver.initialize();

    gsInfo << "Running the transient simulation...\n";
    clock.restart();
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));
        timeSolver.makeTimeStep(timeStep);
        assembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"flappingBeam_CFD3",collection,i + 1 + index_t(warmUpTimeSpan/warmUpTimeStep),numPlotPoints);
        if (validate)
            validation(file,assembler,warmUpTimeSpan + timeStep*(i+1),velocity,pressure);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"flappingBeam_CFD3.pvd\" in Paraview for visualization.\n";
    }
    if (validate)
    {
        file.close();
        gsInfo << "Drag and lift over time are saved to \"flappingBeam_CFD3.txt\".\n";
    }

    return 0;
}
