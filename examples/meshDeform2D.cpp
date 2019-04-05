/// This is an example of generating an isogeometric parametrization using mesh deformation technique
#include <gismo.h>
#include <gsElasticity/gsElMeshing.h>
#include <gsElasticity/gsElasticityFunctions.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Generating isogeometric parametrization by mesh deformation in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/puzzle3_bdry.xml";
    std::string filenameInit = "";
    index_t numSteps = 5;
    real_t poissRatio = 0.45;
    index_t fittingDegree = 0;
    index_t numAdditionalPoints = 0;
    index_t numUniRef = 0;
    index_t numDegreeElev = 0;
    index_t numPlotPoints = 0;
    real_t quality = 0.5;
    index_t maxAdapt = 10;
    index_t maxNewtonIter = 50;
    bool nonLin = true;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Generating isogeometric parametrization by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addString("d","domain","File with an initial domain, if exists",filenameInit);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    cmd.addInt("i","iter","Number of incremental steps for deformation",numSteps);
    cmd.addInt("f","fitDeg","Degree of the fitting curve for simplification",fittingDegree);
    cmd.addInt("a","acc","Number of control points above minimum for curve simplification",numAdditionalPoints);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("e","elev","Number of degree elevetation application",numDegreeElev);
    cmd.addInt("s","sample","Number of points to plot the Jacobain determinant (don't plot if 0)",numPlotPoints);
    cmd.addReal("q","quality","Quality threshold for adaptive incremental loading",quality);
    cmd.addInt("x","maxadapt","Max number of adaptive stepsize halving",maxAdapt);
    cmd.addInt("m","maxiter","Max number of Newton's iterations",maxNewtonIter);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // a set of 4 compatible boundary curves ordered "west-east-south-north"
    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);

    // file name for output
    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\"));
    filename = filename.substr(0,filename.find_last_of("_\\"));

    //=====================================//
                // Algorithm //
    //=====================================//

    for (index_t i = 0; i < numDegreeElev; ++i)
        bdry.degreeElevate();

    for (index_t i = 0; i < numUniRef; ++i)
        bdry.uniformRefine();

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
    gsInfo << "The initial domain is saved to \"" << filename << "_2D_init.xml\".\n";
    gsWrite(initGeo,filename + "_2D_init");

    gsInfo << "Initialized a 2D problem with " << initGeo.patch(0).coefsSize() * 2 << " dofs.\n";

    // boundary condition info
    gsBoundaryConditions<> bdryCurves;
    for (index_t p = 0; p < bdry.nPatches(); ++p)
        bdryCurves.addCondition(0,p+1,condition_type::dirichlet,&(bdry.patch(p)));


    std::vector<gsMultiPatch<> > displacements;
    gsInfo << "Mesh deformation...\n";
    int res = computeMeshDeformation(displacements,initGeo,bdryCurves,poissRatio,numSteps,maxAdapt,quality,nonLin,1e-12,maxNewtonIter);
    (void)res;
    gsInfo << "Plotting the result to the Paraview file \"" << filename << "_2D.pvd\"...\n";
    plotDeformation(displacements,initGeo,filename + "_2D",numPlotPoints);

    // construct deformed geometry
    gsMultiPatch<> geo;
    geo.addPatch(initGeo.patch(0).clone());
    geo.patch(0).coefs() += displacements.back().patch(0).coefs();
    geo.computeTopology();
    gsInfo << "The result of the incremental algorithm is saved to \"" << filename << "_2D.xml\".\n";
    gsWrite(geo,filename + "_2D");

    return 0;
}
