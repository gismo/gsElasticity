/// This is an example of using the nonlinear elasticity solver on a 3D single-patch geometry
///
/// Authors: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Testing the nonlinear elasticity solver in 3D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"spring.xml";
    real_t youngsModulus = 1.0e6;
    real_t poissonsRatio = 0.3;
    index_t materialLaw = material_law::neo_hooke_ln;
    real_t pullZ = 2.0;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numUniRefX = 3;
    index_t numPlotPoints = 10000;
    bool plotMesh = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the linear elasticity solver in 3D.");
    cmd.addReal("z","zpull","Vertical displacement of the top end",pullZ);
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - neoHookeLn, 2 - neoHookeQuad",materialLaw);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("x","xrefine","Number of uniform refinement application in the helical direction",numUniRefX);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addSwitch("m","mesh","Plot computational mesh",plotMesh);
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
    for (index_t i = 0; i < numUniRefX; ++i)
        static_cast<gsTensorNurbsBasis<3,real_t> &>(basis.basis(0)).knots(0).uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> g(0.,0.,0.,3);
    // top end displacement
    gsConstantFunction<> endDisp(pullZ,3);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BC are imposed separately for every component (coordinate)
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,0,2);
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,0);
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,0,1);
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,&endDisp,2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
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
    // constructing stress tensor
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solutionNonlinear,stresses,stress_type::von_mises);

    if (numPlotPoints > 0)
    {
        gsField<> nonlinearSolutionField(geometry,solutionNonlinear);
        gsField<> linearSolutionField(geometry,solutionLinear);
        gsField<> stressField(assembler.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation (nonlinElast)"] = &nonlinearSolutionField;
        fields["Deformation (linElast)"] = &linearSolutionField;
        fields["VonMises"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"spring",numPlotPoints,plotMesh);
        gsInfo << "Open \"spring.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
