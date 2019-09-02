/// This is an example of using the time-dependent Navier-Stokes solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <fstream>

using namespace gismo;

void validation(std::ofstream & ofs,
                const gsNsAssembler<real_t> & assembler, real_t maxInflow,
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
    real_t U_mean = maxInflow * 2./3.; // mean velocity
    real_t drag =  2.*force.at(0)/L/pow(U_mean,2);
    real_t lift =  2.*force.at(1)/L/pow(U_mean,2);

    // evaluating pressure difference at the far front and the far rear points of the cylinder
    gsMatrix<> point(2,1);
    point << 0.5, 0;
    real_t pressureDiff =  pressure.patch(0).eval(point)(0,0) - pressure.patch(2).eval(point)(0,0);

    ofs << drag << " " << lift << " " << pressureDiff << std::endl;
}


int main(int argc, char* argv[]){
    gsInfo << "Testing the time-dependent Navier-Stokes solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 1000;
    real_t viscosity = 0.001;
    real_t density = 1.;
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    real_t maxInflow = 1.5;
    bool subgrid = false;
    bool supg = false;
    real_t timeStep = 0.01;
    real_t timeSpan = 5.;
    bool validate = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the time-dependent Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("d","density","Density of the fluid",density);

    cmd.addSwitch("x","validate","Validate the solver by measuring the drag/lift/pressure difference",validate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
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
    gsMatrix<> box(2,2);
    box << 0.,0.,0.,0.2;
    for (index_t i = 0; i < numBLRef; ++i)
        for (index_t p = 0; p < 4; ++p)
        {
            basisVelocity.refine(p,box);
            basisPressure.refine(p,box);
        }
    // additional velocity refinement for stable mixed FEM
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    gsNsAssembler<real_t> stiffAssembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    stiffAssembler.options().setReal("Viscosity",viscosity);
    stiffAssembler.options().setInt("DirichletValues",dirichlet::interpolation);
    stiffAssembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized system with " << stiffAssembler.numDofs() << " dofs.\n";

    gsMassAssembler<real_t> massAssembler(geometry,basisVelocity,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    gsNsTimeIntegrator<real_t> timeSolver(stiffAssembler,massAssembler);
    timeSolver.setInitialSolution(gsMatrix<>::Zero(stiffAssembler.numDofs(),1));
    timeSolver.initialize();

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
    gsParaviewCollection collection("NS_aroundCylinder");
    // plotting initial condition
    stiffAssembler.constructSolution(gsMatrix<>::Zero(stiffAssembler.numDofs(),1),velocity,pressure);
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,0,numPlotPoints);
    // output file for solver validation
    std::ofstream ofs;
    ofs.open ("validation.txt");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
          // Solving for initial conditions //
    //=============================================//


    timeSolver.options().setInt("Scheme",time_integration_NS::implicit_euler_oseen);
    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    for (index_t i = 0; i < 30; ++i)
    {
        bar.display(i+1,30);
        timeSolver.makeTimeStep(0.1);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                  // Main simulation //
    //=============================================//


    timeSolver.options().setInt("Scheme",time_integration_NS::implicit_euler_oseen);
    clock.restart();
    gsInfo << "Running the simulation with a fine time step...\n";
    for (index_t i = 0; i < timeSpan/timeStep; ++i)
    {
        bar.display(i+1,timeSpan/timeStep);
        timeSolver.makeTimeStep(timeStep);
        stiffAssembler.constructSolution(timeSolver.solutionVector(),velocity,pressure);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"NS_aroundCylinder",collection,i+11,numPlotPoints);
        if (validate)
            validation(ofs,stiffAssembler,maxInflow,velocity,pressure);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
                // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "The results are saved to the Paraview file \"NS_aroundCylinder.pvd\"\n";
    }
    ofs.close();

    return 0;
}
