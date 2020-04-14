/// This is an example of using the nonlinear elasticity solver on a 3D multi-patch geometry
/// The problems is part of the EU project "Terrific".
///
/// Authors: O. Weeger (2012-1015, TU Kaiserslautern),
///          A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"terrific.xml";
    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;
    index_t materialLaw = material_law::saint_venant_kirchhoff;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - NeoHooke_ln",materialLaw);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
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

    // source function, rhs
    gsConstantFunction<> f(0.,0.,0.,3);
    // surface load, neumann BC
    gsConstantFunction<> g(15e7, -10.5e7, 0,3);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BC are imposed separately for every component (coordinate)
    for (index_t d = 0; d < 3; d++)
    {
        bcInfo.addCondition(0,boundary::back,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::back,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
    }
    // Neumann BC are imposed as one function
    bcInfo.addCondition(13,boundary::front,condition_type::neumann,&g);
    bcInfo.addCondition(14,boundary::north,condition_type::neumann,&g);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,f);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("MaterialLaw",materialLaw);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative<real_t> newton(assembler);
    newton.options().setInt("MaxIters",50);
    newton.options().setReal("AbsTol",1e-12);
    newton.options().setInt("Verbosity",solver_verbosity::all);
    newton.options().setInt("Solver",linear_solver::LDLT);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();

    // make the first iteration by hand to get the linear solution
    newton.compute();
    gsInfo << newton.status() << std::endl;
    gsVector<> linSolVector = newton.solution();
    // continue iterations till convergence
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                  // Output //
    //=============================================//

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<> solutionLinear;
    assembler.constructSolution(linSolVector,newton.allFixedDofs(),solutionLinear);
    gsMultiPatch<> solutionNonlinear;
    assembler.constructSolution(newton.solution(),newton.allFixedDofs(),solutionNonlinear);

    if (numPlotPoints > 0)
    {
        gsField<> nonlinearSolutionField(geometry,solutionNonlinear);
        gsField<> linearSolutionField(geometry,solutionLinear);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation (nonlinElast)"] = &nonlinearSolutionField;
        fields["Deformation (linElast)"] = &linearSolutionField;
        gsWriteParaviewMultiPhysics(fields,"terrific",numPlotPoints);
        gsInfo << "Open \"terrific.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
