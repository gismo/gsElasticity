/// This is an example of using the nonlinear elasticity solver on a 3D multi-patch geometry
/// The problems is part of the EU project "Terrific".
///
/// Authors: O. Weeger (2012-1015, TU Kaiserslautern),
///          A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsCompositeMaterial.h>
#include <gsElasticity/gsCompositeMatrix.cpp>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    real_t youngsModulus = 74e9;
    real_t poissonsRatio = 0.33;
    index_t materialLaw = material_law::saint_venant_kirchhoff;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - neoHookeLn, 2 - neoHookeQuad",materialLaw);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> patch, geometry;
    patch.addPatch(gsNurbsCreator<>::BSplineCube());
    patch.patch(0).coefs().col(0)*=10;
    geometry = patch.uniformSplit(0);
    geometry.computeTopology();
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
    gsConstantFunction<> g(0, 0, -1e8,3);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BC are imposed separately for every component (coordinate)
    for (index_t d = 0; d < 3; d++)
    {
        bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
    }
    // Neumann BC are imposed as one function
    bcInfo.addCondition(1,boundary::east,condition_type::neumann,&g);

    //=============================================//
                  // Solving //
    //=============================================//

    gsLinearMaterial<real_t> materialMat(youngsModulus,poissonsRatio,3);

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,f,&materialMat);
    // gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,f);
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
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                  // Output //
    //=============================================//

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<> solutionNonlinear;
    assembler.constructSolution(newton.solution(),newton.allFixedDofs(),solutionNonlinear);
    // constructing stresses
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solutionNonlinear,stresses,stress_components::von_mises);

    if (numPlotPoints > 0)
    {
        gsField<> nonlinearSolutionField(geometry,solutionNonlinear);
        gsField<> stressField(assembler.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &nonlinearSolutionField;
        fields["von Mises"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"terrific_nle",numPlotPoints);
        gsInfo << "Open \"terrific_nle.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
