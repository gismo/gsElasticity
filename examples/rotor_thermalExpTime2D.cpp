/// This is an example of using the thermal expansion solver in a time-dependent setting on a 2D geometry.
/// The problems is part of the EU project "MOTOR".
/// The heat equation is solved using the Crank-Nicolson method.
/// At each time step, the current temperature distribution is used as an input for
/// the stationary thermo-elasticity equation to compute the thermal expansion of the body.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsThermoAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsInfo << "Testing the thermo-elasticity solver in a time-dependent setting in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/rotor_2D.xml";
    real_t fluxValue = 100.; // heat flux on the north boundary
    real_t thExpCoef = 2e-4; // thermal expansion coeffcient of the material
    real_t initTemp = 20.; // initial temperature
    // spatial discretization
    index_t numUniRef = 1;
    index_t numDegElev = 0;
    // time integration
    real_t timeSpan = 1.;
    real_t timeStep = 0.01;
    // output
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the thermo-elasticity solver in a time-dependent setting in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("t","time","Lenght of time integration period",timeSpan);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
         // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basis(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basis.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // heat source function, rhs for the heat equation
    gsConstantFunction<> heatSource(0.,2);
    // boundary temperature, dirichlet BC for the heat equation
    gsConstantFunction<> bTemp(initTemp,2);
    // boundary flux, nuemann BC for the heat equation
    gsConstantFunction<> heatFlux(fluxValue,2);
    // boundary conditions for the heat equation
    gsBoundaryConditions<> bcTemp;
    bcTemp.addCondition(0,boundary::south,condition_type::dirichlet,&bTemp);
    bcTemp.addCondition(0,boundary::north,condition_type::neumann,&heatFlux);

    // gravity, rhs for the linear elasticity equation
    gsConstantFunction<> gravity(0.,0.,2);
    // boundary conditions for the linear elasticity equation
    gsBoundaryConditions<> bcElast;
    // Dirichlet BC are imposed separately for every component (coordinate)
    for (index_t d = 0; d < 2; ++d)
    {   // 0 refers to a patch number, nullptr means that the dirichlet BC is homogeneous, d stats for the displacement component
        bcElast.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,d);
        bcElast.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
        bcElast.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,d);
    }

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // solution fields
    gsMultiPatch<> temperature, displacement;

    // creating an assembler for the heat equation
    gsPoissonAssembler<> heatAssembler(geometry,basis,bcTemp,heatSource);
    gsHeatEquation<real_t> heatSolver(heatAssembler);
    heatSolver.setTheta(0.5);
    gsInfo << "Initialized heat equation system with " << heatAssembler.numDofs() << " dofs.\n";

    // creating elasticity assembler
    gsThermoAssembler<real_t> elastAssembler(geometry,basis,bcElast,gravity,temperature);
    elastAssembler.options().setReal("InitTemp",initTemp);
    elastAssembler.options().setReal("ThExpCoef",thExpCoef);
    gsInfo << "Initialized thermal expansion system with " << elastAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> tempField(geometry,temperature);
    gsField<> dispField(geometry,displacement);

    // setting up Paraview output
    std::map<std::string,const gsField<> *> fields;
    // plotting initial conditions to Paraview
    fields["Temperature"] = &tempField;
    fields["Displacement"] = &dispField;
    gsParaviewCollection collection("rotor");

    gsStopwatch totalClock, iterClock;
    gsProgressBar bar;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    iterClock.restart();
    gsInfo<<"Assembling...\n";
    heatSolver.assemble();
    gsInfo << "Assembled the heat equation system with "
           << heatAssembler.numDofs() << " dofs in " << totalClock.stop() << "s.\n";
    // setting initial conditions
    gsVector<> solVectorTemp;
    solVectorTemp.setOnes(heatAssembler.numDofs());
    solVectorTemp *= initTemp;
    // constructing solution as an IGA function
    heatAssembler.constructSolution(solVectorTemp,temperature);


    totalClock.restart();
    gsInfo<< "Assembling...\n";
    elastAssembler.assemble();
    gsInfo << "Assembled the thermal expansion system with "
           << elastAssembler.numDofs() << " dofs in " << totalClock.stop() << "s.\n";
    // setting initial conditions
    gsVector<> solVectorElast;
    solVectorElast.setZero(elastAssembler.numDofs());
    // constructing solution as an IGA function
    elastAssembler.constructSolution(solVectorElast,displacement);

    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"rotor",collection,0,numPlotPoints);

    //=============================================//
                  // Solving //
    //=============================================//

    real_t timeTemp = 0.;
    real_t timeElast = 0.;

    totalClock.restart();
    gsInfo << "Running the simulation...\n";
    for (int i = 0; i < (index_t)(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,(index_t)(timeSpan/timeStep));

        iterClock.restart();
        // prepairing the matrix for the next time step
        heatSolver.nextTimeStep(solVectorTemp, timeStep);
#ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLDLT solverHeat(heatSolver.matrix());
        solVectorTemp = solverHeat.solve(heatSolver.rhs());
#else
        gsSparseSolver<>::SimplicialLDLT solverHeat(heatSolver.matrix());
        solVectorTemp = solverHeat.solve(heatSolver.rhs());
#endif
        heatAssembler.constructSolution(solVectorTemp,temperature);
        timeTemp += iterClock.stop();

        iterClock.restart();
        // assembling the thermal contribution to the RHS of the thermal expansion system
        elastAssembler.assembleThermo();
#ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLDLT solverElast(elastAssembler.matrix());
        solVectorElast = solverElast.solve(elastAssembler.rhs());
#else
        gsSparseSolver<>::SimplicialLDLT solverElast(elastAssembler.matrix());
        solVectorElast = solverElast.solve(elastAssembler.rhs());
#endif
        elastAssembler.constructSolution(solVectorElast,displacement);
        timeElast += iterClock.stop();

        // output
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"rotor",collection,i+1,numPlotPoints);
    }

    //=============================================//
                // Final touches //
    //=============================================//

    gsInfo << "Complete in: " << secToHMS(totalClock.stop())
           << ", temperature time: " << secToHMS(timeTemp)
           << ", thermal expansion time: " << secToHMS(timeElast) << std::endl;

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"rotor.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
