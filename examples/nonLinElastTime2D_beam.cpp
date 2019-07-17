/// This is an example of using the time-dependent nonlinear elasticity solver on a 2D geometry
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the time-dependent nonlinear elasticity solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/beam.xml";
    index_t numUniRef = 2; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numPlotPoints = 10000;
    index_t numTimeSteps = 100;
    real_t youngsModulus = 200.;//74e9;
    real_t density = 1.;
    real_t force = 0.5;
    // roughly, period = 2 * pi * sqrt(density)/length/sqrt(YoungsModulus) = 2.22
    real_t timeSpan = 6.66;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("k","krefine","Number of degree elevation application",numKRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addInt("i","iter","Number of time steps",numTimeSteps);
    cmd.addReal("y","young","Young's modulus of the material",youngsModulus);
    cmd.addReal("d","density","Density of the material",density);
    cmd.addReal("f","force","Deflecting force for the initial displacement",force);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,force,2);
    gsConstantFunction<> g0(0.,0.,2);

    // material parameters
    real_t poissonsRatio = 0.33;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);

    //=============================================//
                  // Assembly //
    //=============================================//

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
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();

    //=============================================//
         // Solving for initial configuration //
    //=============================================//

    gsElasticityAssembler<real_t> stiffAssemblerInit(geometry,basis,bcInfo,g);
    stiffAssemblerInit.options().setReal("YoungsModulus",youngsModulus);
    stiffAssemblerInit.options().setReal("PoissonsRatio",poissonsRatio);
    stiffAssemblerInit.options().setInt("DirichletValues",dirichlet::interpolation);

    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    stiffAssemblerInit.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << stiffAssemblerInit.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving for the initial configuration...\n";
    clock.restart();
    gsSparseSolver<>::SimplicialLDLT solver(stiffAssemblerInit.matrix());
    gsVector<> solVector = solver.solve(stiffAssemblerInit.rhs());
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    gsMultiPatch<> displacement;
    stiffAssemblerInit.constructSolution(solVector,displacement);

    gsField<> dispField(geometry,displacement);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;

    gsParaviewCollection collection("beam");
    gsWriteParaviewMultiPhysicsTimeStep(fields,"beam",collection,0,numPlotPoints);

    //=============================================//
                  // Transient simulation //
    //=============================================//

    gsElasticityAssembler<real_t> stiffAssembler(geometry,basis,bcInfo,g0);
    stiffAssembler.options().setReal("YoungsModulus",youngsModulus);
    stiffAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    stiffAssembler.options().setInt("DirichletValues",dirichlet::interpolation);

    gsMassAssembler<real_t> massAssembler(geometry,basis,bcInfo,g);
    massAssembler.options().setReal("Density",density);

    gsElTimeIntegrator<real_t> timeSolver(stiffAssembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    // set initial conditions
    timeSolver.setInitialDisplacement(solVector);
    timeSolver.setInitialVeclocity(gsMatrix<>::Zero(stiffAssembler.numDofs(),1));
    timeSolver.initialize();

    gsProgressBar bar;
    real_t timeStep = timeSpan/numTimeSteps;
    clock.restart();
    for (index_t i = 0; i < numTimeSteps; ++i)
    {
        bar.display(i+1,numTimeSteps);
        timeSolver.makeTimeStep(timeStep);
        stiffAssembler.constructSolution(timeSolver.displacementVector(),displacement);
        gsWriteParaviewMultiPhysicsTimeStep(fields,"beam",collection,i+1,numPlotPoints);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";
    collection.save();
    gsInfo << "The results are saved to the Paraview file \"beam.pvd\"\n";
    return 0;
}
