/// This is the 2D linear elasticity benchmark "Infinite plate with circular hole"
/// as described in V.P.Nguyen, C.Anitescu, S.P.A.Bordas, T.Rabczuk, 2015
/// "Isogeometric analysis: An overview and computer implementation aspects".
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is the 2D linear elasticity benchmark: infinite plate with circular hole.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/plateWithHole.xml";
    index_t numUniRef = 5;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;
    bool plotMesh = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is the 2D linear elasticity benchmark: infinite plate with circular hole.");
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
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

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsFunctionExpr<> analyticalStresses("1-1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))",2);
    // boundary load neumann BC
    gsFunctionExpr<> traction("(-1+1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (y==4)",
                              "(1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (y==4)",2);
    // material parameters
    real_t youngsModulus = 1.0e3;
    real_t poissonsRatio = 0.3;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::neumann,&traction);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    gsLinearMaterial<real_t> materialMat(youngsModulus,poissonsRatio);

    // creating assembler
    // gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);//,materialMat);
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g,&materialMat);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
#endif

    //=============================================//
                      // Output //
    //=============================================//

    // constructing displacement as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),solution);
    // constructing stress tensor
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_components::all_2D_vector);

    if (numPlotPoints > 0)
    {
        // constructing an IGA field (geometry + solution) for displacement
        gsField<> solutionField(assembler.patches(),solution);
        // constructing an IGA field (geometry + solution) for stresses
        gsField<> stressField(assembler.patches(),stresses,true);
        // analytical stresses
        gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Stress"] = &stressField;
        fields["StressAnalytical"] = &analyticalStressField;
        gsWriteParaviewMultiPhysics(fields,"plateWithHole",numPlotPoints,plotMesh);
        gsInfo << "Open \"plateWithHole.pvd\" in Paraview for visualization. Stress wiggles on the left side are caused by "
                  "a singularity in the parametrization.\n";
    }

    // eval stress at the top of the circular cut
    gsMatrix<> A(2,1);
    A << 1.,0.; // parametric coordinates for the isogeometric solution
    gsMatrix<> res;
    stresses.piece(0).eval_into(A,res);
    A << 0., 1.; // spatial coordinates for the analytical solution
    gsMatrix<> analytical;
    analyticalStresses.eval_into(A,analytical);
    gsInfo << "XX-stress at the top of the circle: " << res.at(0) << " (computed), " << analytical.at(0) << " (analytical)\n";
    gsInfo << "YY-stress at the top of the circle: " << res.at(1) << " (computed), " << analytical.at(1) << " (analytical)\n";
    gsInfo << "XY-stress at the top of the circle: " << res.at(2) << " (computed), " << analytical.at(2) << " (analytical)\n";


    return 0;
}
