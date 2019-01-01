/// This is an example of generating an isogeometric parametrization using mesh deformation technique
#include <gismo.h>
#include <gsElasticity/gsElMeshing.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "Isogeometric mesh generation by mesh deformation in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/puzzle4_bdry.xml";
    std::string filenameInit = "";
    int numSteps = 3;
    real_t poissRatio = 0.45;
    int numAdditionalPoints = 0;
    int numUniRef = 0;
    int numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Isogeometric mesh generation by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addString("d","domain","File with an initial domain, if exists",filenameInit);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    cmd.addInt("i","ites","Number of incremental steps for deformation",numSteps);
    cmd.addInt("a","acc","Number of control points above minimum for curve simplification",numAdditionalPoints);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // a set of 4 compatible boundary curves ordered "west-east-south-north"
    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);

    //=====================================//
                // Algorithm //
    //=====================================//

    for (int i = 0; i < numUniRef; ++i)
        bdry.uniformRefine();

    // simplifying the boundary curves
    gsMultiPatch<> simpleBdry;
    for (int p = 0; p < bdry.nPatches(); ++p)
        simpleBdry.addPatch(simplifyCurve(bdry.patch(p),numAdditionalPoints,1000));

    // creating a coons patch out of simplified boundary curves to serve as an initial domain
    gsCoonsPatch<real_t> coonsPatch(simpleBdry);
    coonsPatch.compute();
    gsMultiPatch<> initGeo;
    initGeo.addPatch(coonsPatch.result());
    initGeo.computeTopology();

    // boundary condition info
    gsBoundaryConditions<> bcDeform;
    for (int p = 0; p < bdry.nPatches(); ++p)
        bcDeform.addCondition(0,p+1,condition_type::dirichlet,&(bdry.patch(p)));

    std::vector<std::vector<gsMatrix<> > > result;
    computeDeformationInc(result,initGeo,bcDeform,numSteps,poissRatio);

    //=====================================//
                // Output //
    //=====================================//

    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\"));
    filename = filename.substr(0,filename.find_last_of("_\\"));

    plotDeformation(result,initGeo,filename + "_2D",true,numPlotPoints);

    gsMultiPatch<> geo;
    applyDeformation(result,initGeo,geo);
    gsWrite(geo,filename + "_2D");

    return 0;
}
