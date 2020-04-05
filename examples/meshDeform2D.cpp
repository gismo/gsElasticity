/// This is an example of generating an isogeometric parametrization using mesh deformation technique
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsALE.h>
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
    real_t poissRatio = 0.45;
    real_t stiffDegree = 0.;
    index_t ALEmethod = ale_method::TINE;
    /// Output options
    index_t numPlotPoints = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Generating isogeometric parametrization by mesh deformation in 2D.");
    cmd.addPlainString("name","File with boundary curves.",filename);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    /// Initial domain options
    cmd.addInt("f","fitDeg","Polynomial degree of the coarse fitting curve",fittingDegree);
    cmd.addInt("c","acp","Number of additional control points above minimal for the coarse fitting curve",numAdditionalPoints);
    cmd.addString("d","domain","Optional input file with a before-hand prepaired initial domain",filenameInit);
    /// Deformation options
    cmd.addInt("i","iter","Number of incremental loading steps used during deformation",numSteps);
    cmd.addReal("p","poiss","Poisson's ratio for the elasticity model",poissRatio);
    cmd.addReal("x","xjac","Stiffening degree for the Jacobian-based local stiffening",stiffDegree);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
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
    initGeo.computeTopology();

    gsInfo << "Initialized a 2D problem with " << initGeo.patch(0).coefsSize() * 2 << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // container for displacements
    std::vector<gsMultiPatch<> > displacements;
    gsStopwatch clock;
    gsProgressBar bar;

    //=============================================//
             // Setting mesh deformation //
    //=============================================//

    // MP displacement object that defines bdry displacement; interior can be arbitrary (e.g., Coons patch)
    gsMultiPatch<> bdryDisplacement;
    gsCoonsPatch<real_t> coonsPatch(bdry);
    coonsPatch.compute();
    bdryDisplacement.addPatch(coonsPatch.result());
    bdryDisplacement.patch(0).coefs() -= initGeo.patch(0).coefs();
    bdryDisplacement.patch(0).coefs() /= numSteps;
    // Boundary sides to deform
    gsInterfaceFSI interface;
    interface.addSide(0,boundary::south,0,boundary::south);
    interface.addSide(0,boundary::north,0,boundary::north);
    interface.addSide(0,boundary::west,0,boundary::west);
    interface.addSide(0,boundary::east,0,boundary::east);
    // mesh deformation object
    gsALE<real_t> meshDeformer(initGeo,bdryDisplacement,interface,ale_method::method(ALEmethod));
    meshDeformer.options().setReal("LocalStiff",stiffDegree);
    meshDeformer.options().setReal("PoissonsRatio",poissRatio);

    //=====================================//
                  // Deforming //
    //=====================================//

    gsInfo << "Solving...\n";
    clock.restart();
    for (index_t i = 0; i < numSteps; ++i)
    {
        bar.display(i+1,numSteps);
        // deform mesh to match the current bdry displacement
        meshDeformer.updateMesh();
        // save the displacement
        displacements.push_back(gsMultiPatch<>());
        meshDeformer.constructSolution(displacements.back());
        // increase the bdry displacement for the next step
        bdryDisplacement.patch(0).coefs() *= 1.*(i+2)/(i+1);
    }

    gsInfo << "Solved in "<< clock.stop() <<"s.\n";

    //=====================================//
                // Output //
    //=====================================//

    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\")); // file name without an extension
    filename = filename.substr(0,filename.find_last_of("_\\"));

    // save initial domain
    gsInfo << "The initial domain is saved to \"" << filename << "_2D_init.xml\".\n";
    gsWrite(initGeo,filename + "_2D_init");
    // save resulting domain
    gsInfo << "The result of the deformation algorithm is saved to \"" << filename << "_2D.xml\".\n";
    gsWrite(initGeo,filename + "_2D");
    // plotting all intermediate deformed meshes
    plotDeformation(initGeo,displacements,filename,numPlotPoints);
    gsInfo << "Open \"" << filename << "_mesh.pvd\" in Paraview for visualization.\n";

    return 0;
}
