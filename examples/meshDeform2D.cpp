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

    std::string filename = ELAST_DATA_DIR"/rotor_bdry.xml";
    std::string filenameInit = "";
    index_t numSteps = 3;
    real_t poissRatio = 0.45;
    index_t fittingDegree = 0;
    index_t numAdditionalPoints = 0;
    index_t materialLaw = 1;
    index_t numUniRef = 0;
    index_t numDegreeElev = 0;
    index_t numPlotPoints = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Generating isogeometric parametrization by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addString("d","domain","File with an initial domain, if exists",filenameInit);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    cmd.addInt("i","ites","Number of incremental steps for deformation",numSteps);
    cmd.addInt("f","fitDeg","Degree of the fitting curve for simplification",fittingDegree);
    cmd.addInt("a","acc","Number of control points above minimum for curve simplification",numAdditionalPoints);
    cmd.addInt("l","law","Material law: 0 - St.V.-K., 1 - NeoHooke_ln, 2 - NeoHooke_2; if not set, no nonlin solution",materialLaw);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("e","elev","Number of degree elevetation application",numDegreeElev);
    cmd.addInt("s","sample","Number of points to plot the Jacobain determinant (don't plot if 0)",numPlotPoints);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // a set of 4 compatible boundary curves ordered "west-east-south-north"
    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);

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
    gsInfo << "Initialized a 2D problem with " << initGeo.patch(0).coefsSize() * 2 << " dofs.\n";

    // boundary condition info
    gsBoundaryConditions<> bdryCurves;
    for (index_t p = 0; p < bdry.nPatches(); ++p)
        bdryCurves.addCondition(0,p+1,condition_type::dirichlet,&(bdry.patch(p)));

    gsInfo << "Computing deformation using linear elasticity with incremental Dirichelt BC...\n";
    std::vector<std::vector<gsMatrix<> > > deformation;
    computeDeformationPlus(deformation,initGeo,bdryCurves,numSteps,poissRatio);

    // construct deformed geometry
    gsMultiPatch<> geo;
    applyDeformation(deformation,initGeo,geo);

    // construct displacement field
    gsMultiPatch<> displacement;
    constructDisplacement(deformation,initGeo,displacement);

    gsMultiPatch<> nonlinGeo;
    if (materialLaw >= 0)
    {
        gsInfo << "Solving a nonlinear problem using the incremental solution as an initial guess...\n";
        computeDeformationNonlin(nonlinGeo,initGeo,bdryCurves,displacement,materialLaw,poissRatio);
    }

    //=====================================//
                // Output //
    //=====================================//

    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\"));
    filename = filename.substr(0,filename.find_last_of("_\\"));

    gsInfo << "The initial domain is saved to \"" << filename << "_2D_init.xml\".\n";
    gsWrite(initGeo,filename + "_2D_init");

    gsInfo << "Plotting the result of the incremental algorithm to the Paraview file \"" << filename << "_2D_lin.pvd\"...\n";
    plotDeformation(deformation,initGeo,filename + "_2D_lin",numPlotPoints);
    gsInfo << "The result of the incremental algorithm is saved to \"" << filename << "_2D_lin.xml\".\n";
    gsWrite(geo,filename + "_2D_lin");

    if (materialLaw >= 0)
    {
        gsInfo << "Plotting the result of the nonlinear algorithm to the Paraview file \"" << filename << "_2D_nl.pvd\"...\n";
        plotGeometry(nonlinGeo,filename + "_2D_nl",numPlotPoints);
        gsInfo << "The result of the nonlinear algorithm is saved to \"" << filename << "_2D_nl.xml\".\n";
        gsWrite(nonlinGeo,filename + "_2D_nl");       
    }

    return 0;
}
