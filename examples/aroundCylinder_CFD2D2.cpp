/// This is the incompressible Navier-Stokes solver benchmark 2D-2 from this project:
/// http://www.featflow.de/en/benchmarks/cfdbenchmarking.htmly
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <fstream>

using namespace gismo;

void refineBoundaryLayer(gsMultiBasis<> & velocity, gsMultiBasis<> & pressure)
{
    gsMatrix<> boxSouth(2,2);
    boxSouth << 0.,0.,0.,0.2;
    for (index_t p = 0; p < 4; ++p)
    {
        velocity.refine(p,boxSouth);
        pressure.refine(p,boxSouth);
    }
}

void validation(std::ofstream & ofs,
                const gsNsAssembler<real_t> & assembler, real_t maxInflow, real_t time,
                const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure)
{
    // computing force acting on the surface of the cylinder
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);
    real_t L = 0.1; // characteristic lenght
    real_t drag =  2.*force.at(0)/L/pow(maxInflow * 2./3.,2);
    real_t lift =  2.*force.at(1)/L/pow(maxInflow * 2./3.,2);

    // evaluating pressure difference at the far front and the far rear points of the cylinder
    gsMatrix<> point(2,1);
    point << 0.5, 0;
    real_t pressureDiff =  pressure.patch(0).eval(point)(0,0) - pressure.patch(2).eval(point)(0,0);
    // print time-drag-lift-presDiff
    ofs << time << " " << drag << " " << lift << " " << pressureDiff << std::endl;
}


int main(int argc, char* argv[]){
    gsInfo << "Benchmark CFD2D2: transient flow of an incompressible fluid.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numPlotPoints = 900;
    real_t viscosity = 0.001;
    real_t density = 1.;
    real_t maxInflow = 1.5;
    bool subgrid = false;
    real_t timeSpan = 5;
    real_t timeStep = 0.01;
    bool validate = true;
    real_t warmUpTimeSpan = 2.;
    real_t warmUpTimeStep = 0.1;
    real_t theta = 0.5;


    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the time-dependent Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addSwitch("x","validate","Validate the solver by measuring the drag/lift/pressure difference",validate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
          // Setting assemblers and solvers //
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
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    // inflow velocity profile U(y) = U_max*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,&inflow,0);
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
    }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // creating assembler
    gsNsAssembler<real_t> stiffAssembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    stiffAssembler.options().setReal("Viscosity",viscosity);
    stiffAssembler.options().setReal("Density",density);
    gsInfo << "Initialized system with " << stiffAssembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsNsTimeIntegrator<real_t> timeSolver(stiffAssembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration_NS::theta_scheme_linear);
    timeSolver.options().setReal("Theta",theta);
    timeSolver.options().setInt("Verbosity",newton_verbosity::all);

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // containers for solution as IGA functions
    gsMultiPatch<> velocity, pressure;
    // isogeometric fields (geometry + solution)
    gsField<> velocityField(geometry,velocity);
    gsField<> pressureField(geometry,pressure);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;
    // paraview collection of time steps
    gsParaviewCollection collection("aroundCylinder_CFD2D2");

    // output file for solver validation
    std::ofstream file;
    if (validate)
        file.open ("aroundCylinder_CFD2D2.txt");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
                   // Warming up //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    stiffAssembler.getFixedDofs(0,boundary::north,inflowDDoFs);

    // set all Dirichlet DoFs to zero
    stiffAssembler.homogenizeFixedDofs(-1);
    gsMatrix<> solVector = gsMatrix<>::Zero(stiffAssembler.numDofs(),1);

    // set initial velocity: zero free and fixed DoFs
    timeSolver.setSolutionVector(solVector);
    timeSolver.setFixedDofs(stiffAssembler.allFixedDofs());
    timeSolver.initialize();

    // consruct and plot initial velocity
    stiffAssembler.constructSolution(solVector,stiffAssembler.allFixedDofs(),velocity,pressure);
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"aroundCylinder_CFD2D2",collection,0,numPlotPoints);

    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    for (index_t i = 0; i < index_t(warmUpTimeSpan/warmUpTimeStep); ++i)
    {
        bar.display(i+1,index_t(warmUpTimeSpan/warmUpTimeStep));
        stiffAssembler.setFixedDofs(0,boundary::north,inflowDDoFs*(1-cos(M_PI*warmUpTimeStep*(i+1)/warmUpTimeSpan))/2);
        timeSolver.makeTimeStep(warmUpTimeStep);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),timeSolver.allFixedDofs(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"aroundCylinder_CFD2D2",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Running the main simulation...\n";
    clock.restart();
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));
        timeSolver.makeTimeStep(timeStep);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),timeSolver.allFixedDofs(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"aroundCylinder_CFD2D2",collection,i + 1 + index_t(warmUpTimeSpan/warmUpTimeStep),numPlotPoints);
        if (validate)
            validation(file,stiffAssembler,maxInflow,warmUpTimeSpan + timeStep*(i+1),velocity,pressure);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"aroundCylinder_CFD2D2.pvd\" in Paraview for visualization.\n";
    }
    if (validate)
    {
        file.close();
        gsInfo << "Drag and lift coefficients and the pressure difference over time are saved to \"aroundCylinder_CFD2D2.txt\".\n";
    }

    return 0;
}
