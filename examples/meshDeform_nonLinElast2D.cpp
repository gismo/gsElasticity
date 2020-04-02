/// This is an example of generating an isogeometric parametrization using mesh deformation technique
/// with incremental deformation based on nonlinear elasticity (TINE method)
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Generating isogeometric parametrization by mesh deformation in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // input file with a set of 4 compatible boundary curves ordered "west-east-south-north"
    std::string filename = ELAST_DATA_DIR"/puzzle3_bdry.xml";
    index_t numUniRef = 0;
    /// Initial domain options
    index_t fittingDegree = 0;
    index_t numAdditionalPoints = 0;
    std::string filenameInit = "";
    /// Deformation options
    index_t numSteps = 5;
    index_t numIter = 1;
    real_t poissRatio = 0.45;
    index_t law = material_law::neo_hooke_ln;
    /// Output options
    index_t numPlotPoints = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Generating isogeometric parametrization by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    /// Initial domain options
    cmd.addInt("f","fitDeg","Polynomial degree of the coarse fitting curve",fittingDegree);
    cmd.addInt("c","acp","Number of additional control points above minimal for the coarse fitting curve",numAdditionalPoints);
    cmd.addString("d","domain","Optional input file with a previously prepaired initial domain",filenameInit);
    /// Deformation options
    cmd.addInt("i","iter","Number of incremental loading steps used during deformation",numSteps);
    cmd.addInt("n","numIter","Number of Newton's iterations at each non-final incremental loading step",numIter);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    cmd.addInt("l","law","Material law used in the material model; 0 - St. Venant-Kirchhoff, 1 - ln-neo-Hooke, 2 - quad-neo-Hooke",law);
    /// Output options
    cmd.addInt("s","sample","Number of points to plot the Jacobain determinant (don't plot if 0)",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=====================================//
          // Initial domain generation //
    //=====================================//

    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);

    for (index_t i = 0; i < numUniRef; ++i)
        bdry.uniformRefine();

    // Initial domain generation
    gsMultiPatch<> initGeo;
    if (filenameInit.empty())
    {
        // simplifying the boundary curves
        gsMultiPatch<> simpleBdry;
        for (size_t p = 0; p < bdry.nPatches(); ++p)
            simpleBdry.addPatch(simplifyCurve(bdry.patch(p),numAdditionalPoints,fittingDegree,1000));

        // creating a coons patch out of simplified boundary curves to serve as an initial domain
        gsCoonsPatch<real_t> coonsPatch(simpleBdry);
        coonsPatch.compute();
        initGeo.addPatch(coonsPatch.result());
    }
    else
        gsReadFile<>(filenameInit,initGeo);

    gsInfo << "Initialized a 2D problem with " << initGeo.patch(0).coefsSize() * 2 << " dofs.\n";

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsMultiBasis<> basis(initGeo);
    // first tell assembler to allocate space for Dirichlet DoFs; we specify the values later
    gsBoundaryConditions<> bcInfo;
    for (index_t s = 1; s < 5; ++s)
    {
        bcInfo.addCondition(0,s,condition_type::dirichlet,0,0);
        bcInfo.addCondition(0,s,condition_type::dirichlet,0,1);
    }
    gsConstantFunction<> g(0.,0.,2);

    //=====================================//
                // Solving //
    //=====================================//

    gsElasticityAssembler<real_t> assembler(initGeo,basis,bcInfo,g);
    assembler.options().setReal("PoissonsRatio",poissRatio);
    assembler.options().setInt("MaterialLaw",law);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating the nonlinear solver
    gsIterative<real_t> newton(assembler);
    newton.options().setInt("Verbosity",solver_verbosity::all);
    newton.options().setInt("MaxIters",numIter);
    newton.options().setInt("Solver",linear_solver::LDLT);

    // container for displacements
    std::vector<gsMultiPatch<> > displacements;

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();

    for (index_t i = 0; i < numSteps; ++i)
    {
        gsInfo << "Loading: " << index_t(100.*i/numSteps) << "% -> " << index_t(100.*(i+1)/numSteps) << "%\n";
        // setting Dirichlet DoFs
        for (index_t s = 1; s < 5; ++s)
            assembler.setFixedDofs(0,s,(bdry.patch(s-1).coefs() - initGeo.patch(0).boundary(s)->coefs())/numSteps);

        if (i == numSteps-1)
            newton.options().setInt("MaxIters",50);
        newton.reset();
        newton.solve();
        displacements.push_back(gsMultiPatch<>());
        assembler.constructSolution(newton.solution(),newton.allFixedDofs(),displacements.back());
    }
    gsInfo << "Solved in "<< clock.stop() <<"s.\n";

    //=====================================//
                // Output //
    //=====================================//

    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\")); // file name without an extension
    filename = filename.substr(0,filename.find_last_of("_\\"));

    // plotting all intermediate deformed meshes
    plotDeformation(initGeo,displacements,filename,numPlotPoints);
    gsInfo << "Open \"" << filename << "_mesh.pvd\" in Paraview for visualization.\n";
    // save initial domain
    gsInfo << "The initial domain is saved to \"" << filename << "_2D_init.xml\".\n";
    gsWrite(initGeo,filename + "_2D_init");
    initGeo.patch(0).coefs() += displacements.back().patch(0).coefs();
    gsInfo << "The result of the deformation algorithm is saved to \"" << filename << "_2D.xml\".\n";
    gsWrite(initGeo,filename + "_2D");

    return 0;
}
