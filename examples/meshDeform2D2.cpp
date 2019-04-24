/// This is an example of generating an isogeometric parametrization using mesh deformation technique

#include <gismo.h>
#include <gsElasticity/gsElMeshing.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityNewton.h>

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
    index_t numDegreeElev = 0;
    /// Initial domain options
    index_t fittingDegree = 0; // polynomial degree of the coarse fitting curve
    index_t numAdditionalPoints = 0; // number of additional control points above minimal for the coarse fitting curve
    std::string filenameInit = ""; // optional input file with a previously prepaired initial domain
    /// Deformation options
    index_t numSteps = 5; // number of incremental loading steps used during deformation
    index_t numIter = 1; // number of Newton's iterations at each non-final incremental loading step
    real_t poissRatio = 0.45;
    /// Bijectivity-preserving adaptivity settings
    index_t maxAdapt = 10; // max number of adaptive stepsize halving applications
    real_t quality = 0.5; // quality ratio to preserve
    /// Output options
    index_t numPlotPoints = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Generating isogeometric parametrization by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("e","elev","Number of degree elevetation application",numDegreeElev);
    /// Initial domain options
    cmd.addInt("f","fitDeg","Polynomial degree of the coarse fitting curve",fittingDegree);
    cmd.addInt("c","acp","Number of additional control points above minimal for the coarse fitting curve",numAdditionalPoints);
    cmd.addString("d","domain","Optional input file with a previously prepaired initial domain",filenameInit);
    /// Deformation options
    cmd.addInt("i","iter","Number of incremental loading steps used during deformation",numSteps);
    cmd.addInt("n","numIter","Number of Newton's iterations at each non-final incremental loading step",numIter);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    /// Bijectivity-preserving adaptivity settings
    cmd.addInt("a","adapt","Max number of adaptive stepsize halving",maxAdapt);
    cmd.addReal("q","quality","Quality threshold for adaptive incremental loading",quality);
    /// Output options
    cmd.addInt("s","sample","Number of points to plot the Jacobain determinant (don't plot if 0)",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=====================================//
                // Algorithm //
    //=====================================//

    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);
    for (index_t i = 0; i < numDegreeElev; ++i)
        bdry.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        bdry.uniformRefine();

    // Initial domain generation
    gsMultiPatch<> initGeo;
    if (filenameInit.empty())
    {
        // simplifying the boundary curves
        gsMultiPatch<> simpleBdry;
        for (index_t p = 0; p < bdry.nPatches(); ++p)
            simpleBdry.addPatch(simplifyCurve(bdry.patch(p),numAdditionalPoints,fittingDegree,1000));

        // creating a coons patch out of simplified boundary curves to serve as an initial domain
        gsCoonsPatch<real_t> coonsPatch(simpleBdry);
        coonsPatch.compute();
        initGeo.addPatch(coonsPatch.result());
    }
    else
        gsReadFile<>(filenameInit,initGeo);

    initGeo.computeTopology();
    gsInfo << "Initialized a 2D problem with " << initGeo.patch(0).coefsSize() * 2 << " dofs.\n";

    gsMultiBasis<> basis(initGeo);
    // first tell assembler to allocate space for Dirichlet DoFs; we specify the values later
    gsBoundaryConditions<> bcInfo;
    for (index_t s = 1; s < 5; ++s)
    {
        bcInfo.addCondition(0,s,condition_type::dirichlet,0,0);
        bcInfo.addCondition(0,s,condition_type::dirichlet,0,1);
    }
    gsConstantFunction<> g(0.,0.,2);

    gsElasticityAssembler<real_t> assembler(initGeo,basis,bcInfo,g);
    assembler.options().setReal("PoissonsRatio",poissRatio);
    assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    for (index_t s = 1; s < 5; ++s)
        assembler.setDirichletDofs(0,s,bdry.patch(s-1).coefs() - initGeo.patch(0).boundary(s)->coefs());

    gsElasticityNewton<real_t> newton(assembler);
    newton.options().setInt("NumIncStep",numSteps);
    newton.options().setInt("MaxIterNotLast",numIter);
    newton.options().setInt("Verbosity",newtonVerbosity::all);
    newton.options().setInt("Save",newtonSave::all);

    newton.solve();

    //=====================================//
                // Output //
    //=====================================//

    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\"));
    filename = filename.substr(0,filename.find_last_of("_\\"));

    gsInfo << "The initial domain is saved to \"" << filename << "_2D_init.xml\".\n";
    gsWrite(initGeo,filename + "_2D_init");

    newton.plotDeformation(initGeo,filename,numPlotPoints);

    initGeo.patch(0).coefs() += newton.solution().patch(0).coefs();
    gsInfo << "The result of the deformation algorithm is saved to \"" << filename << "_2D.xml\".\n";
    gsWrite(initGeo,filename + "_2D");

    return 0;
}
