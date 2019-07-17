/// This is an example of using the thermal expansion solver in a time-dependent setting on a 2D geometry.
/// The heat equation is solved using the Crank-Nicolson method.
/// At time step the current temperature distribution is used as an input for
/// the stationary thermo-elasticity equation to compute the thermal expansion of the body.

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
    int numUniRef = 0; // number of h-refinements
    int numKRef = 0; // number of k-refinements
    int numPlotPoints = 10000;

    real_t endTime = 1.;
    int numSteps = 25;
    real_t theta = 1.;
    real_t fluxValue = 100.;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the thermo-elasticity solver in a time-dependent setting in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addInt("n","num","Number of time steps",numSteps);
    cmd.addReal("e","end","End of the integration period",endTime);
    cmd.addReal("t","theta","Time integration: 0 - explicit Euler, 1 - implicit Euler, 0.5 - Crank-Nicolson",theta);
    cmd.addReal("f","flux","Heat flux at the boundary",fluxValue);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // material parameters
    real_t youngsModulus = 74e9; // doesn't matter for the simulation
    real_t poissonsRatio = 0.33; // doesn't matter for the simulation
    real_t thExpCoef = 2e-4;
    real_t initTemp = 20.;

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
    bcElast.addCondition(0,boundary::south,condition_type::dirichlet,0,0); // the first number is a patch number
    bcElast.addCondition(0,boundary::south,condition_type::dirichlet,0,1); // the second number is a NULL pointer meaning a zero function
    bcElast.addCondition(0,boundary::west,condition_type::dirichlet,0,0);  // the third number is a component (coordinate) number
    bcElast.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    bcElast.addCondition(0,boundary::east,condition_type::dirichlet,0,0);
    bcElast.addCondition(0,boundary::east,condition_type::dirichlet,0,1);

    //=====================================================//
                  // Assembly & output setup //
    //=====================================================//

    gsStopwatch clock;

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating basis
    gsMultiBasis<> basis(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basis.degreeElevate();
        basis.uniformRefine();
    }
    for (int i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    // creating an assembler for the heat equation
    gsPoissonAssembler<> stationary(geometry,basis,bcTemp,heatSource);
    gsHeatEquation<real_t> heatAssembler(stationary);
    heatAssembler.setTheta(theta);

    clock.restart();
    gsInfo<<"Assembling...\n";
    heatAssembler.assemble();
    gsInfo << "Assembled the heat equation system (matrices and a load vector) with "
           << heatAssembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    // setting initial conditions
    gsVector<> solVectorTemp;
    solVectorTemp.setConstant(heatAssembler.numDofs(),initTemp);

    // constructing solution as an IGA function
    gsMultiPatch<> solutionTemp;
    heatAssembler.constructSolution(solVectorTemp,solutionTemp);
    gsField<> tempField(stationary.patches(),solutionTemp);

    // creating elasticity assembler
    gsThermoAssembler<real_t> elastAssembler(geometry,basis,bcElast,gravity,solutionTemp);
    elastAssembler.options().setReal("YoungsModulus",youngsModulus);
    elastAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elastAssembler.options().setReal("InitTemp",initTemp);
    elastAssembler.options().setReal("ThExpCoef",thExpCoef);

    clock.restart();
    gsInfo<< "Assembling...\n";
    elastAssembler.assemble();
    gsInfo << "Assembled the elasticity system (matrix and load vector) with "
           << elastAssembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    // setting initial conditions
    gsMatrix<> solVectorElast;
    solVectorElast.setZero(elastAssembler.numDofs(),1);

    // constructing solution as an IGA function
    gsMultiPatch<> solutionElast;
    elastAssembler.constructSolution(solVectorElast,solutionElast);
    gsField<> elastField(elastAssembler.patches(),solutionElast);

    // setting up Paraview output
    gsParaviewCollection collection("rotor");
    std::map<std::string,const gsField<> *> fields;
    // plotting initial conditions to Paraview
    fields["Temperature"] = &tempField;
    fields["Displacement"] = &elastField;

    gsWriteParaviewMultiPhysicsTimeStep(fields,"rotor",collection,0,numPlotPoints);

    //=====================================================//
                  // Main time loop //
    //=====================================================//

    real_t Dt = endTime / numSteps ;
    gsInfo << "Solving...\n";
    clock.restart();
    gsProgressBar bar;
    for ( int i = 1; i<=numSteps; ++i)
    {
        // display progress bar
        bar.display(i,numSteps);
        // prepairing the matrix for the next time step
        heatAssembler.nextTimeStep(solVectorTemp, Dt);
        // solving the heat system
        gsSparseSolver<>::LU solverHeat(heatAssembler.matrix());
        solVectorTemp = solverHeat.solve(heatAssembler.rhs());
        // constructing solution as an IGA function
        stationary.constructSolution(solVectorTemp,solutionTemp);
        // assembling the thermal contribution to the RHS
        elastAssembler.assembleThermo();
        // solving elasticity system
        gsSparseSolver<>::LU solverElast(elastAssembler.matrix());
        solVectorElast = solverElast.solve(elastAssembler.rhs());
        // constructing solution as an IGA function
        elastAssembler.constructSolution(solVectorElast,solutionElast);
        gsField<> elastField(elastAssembler.patches(),solutionElast);
        // plotting to Paraview
        fields["Temperature"] = &tempField;
        fields["Displacement"] = &elastField;
        gsWriteParaviewMultiPhysicsTimeStep(fields,"rotor",collection,i,numPlotPoints);
    }

    gsInfo << "Total time elapsed: " << clock.stop() << "s\n";
    collection.save();

    return 0;
}
