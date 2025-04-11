/// This is the "Cook's membrane" benchmark solved using the linear elasticity solver with a mixed displacement-pressure formulation.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is Cook's membrane benchmark with mixed linear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = gsElasticity_DATA"/cooks.xml";
    real_t youngsModulus = 240.565e6;
    real_t poissonsRatio = 0.4;
    index_t numUniRef = 3;
    index_t numDegElev = 1;
    bool subgridOrTaylorHood = false;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with mixed linear elasticity solver.");
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation applications",numDegElev);
    cmd.addSwitch("e","element","Mixed element: false = subgrid (default), true = Taylor-Hood",subgridOrTaylorHood);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisDisplacement.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional displacement refinement for stable mixed FEM
    if (!subgridOrTaylorHood) // subgrid
        basisDisplacement.uniformRefine();
    else  // Taylor-Hood
        basisDisplacement.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // neumann BC
    gsConstantFunction<> f(0.,625e4,2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (index_t d = 0; d < 2; ++d)
        bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
    bcInfo.addCondition(0,boundary::east,condition_type::neumann,&f);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler.matrix());
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

    // constructing solution as an IGA function
    gsMultiPatch<> displacement, pressure;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),displacement,pressure);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

    if (numPlotPoints > 0) // visualization
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),displacement);
        gsField<> pressureField(assembler.patches(),pressure);
        gsField<> stressField(assembler.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &displacementField;
        fields["Pressure"] = &pressureField;
        fields["von Mises"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"cooks",numPlotPoints);
        gsInfo << "Open \"cooks.pvd\" in Paraview for visualization.\n";
    }

    // validation
    gsMatrix<> A(2,1);
    A << 1.,1.;
    A = displacement.patch(0).eval(A);
    gsInfo << "X-displacement of the top-right corner: " << A.at(0) << std::endl;
    gsInfo << "Y-displacement of the top-right corner: " << A.at(1) << std::endl;

    return 0;
}
