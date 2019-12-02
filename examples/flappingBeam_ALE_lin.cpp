/// This is the benchmark for computing ALE mapping for an FSI problem based on this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// There is no flow in the benchmark. The flow domain deformation is driven by a freely oscillating beam.
/// The ALE mapping is computed using the linear elasticity method with Jacobian-based local stiffening.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Testing the ALE mapping construction in 2D.\n";

    std::string filenameALE = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t loading = 2.;
    real_t densitySolid = 1.0e3;
    real_t timeStep = 0.01;
    real_t timeSpan = 0.91;
    real_t poissonsRatioMesh = 0.3;
    real_t stiffDegree = 2.3;
    bool stop = true;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the steady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Pois ratio for mesh",poissonsRatioMesh);
    cmd.addReal("l","load","Gravity loading acting on the beam",loading);
    cmd.addReal("x","xjac","Stiffening degree for the Jacobian-based local stiffening",stiffDegree);
    cmd.addSwitch("b","break","Break the simulation if the mesh becomes locally inverted",stop);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);
    gsMultiPatch<> geoFlowFull; // this is a full flow geometry
    gsReadFile<>(filenameALE, geoFlowFull);
    // this is a copy of the flow domain with only the first 6 patches.
    gsMultiPatch<> geoALE;
    for (index_t p = 0; p < 6; ++p)
        geoALE.addPatch(geoFlowFull.patch(p).clone());
    geoALE.computeTopology();

    // creating bases
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        // geoALE is also refined in order to move its control points using a refined diplacement field
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> gBeam(0.,loading*densitySolid,2);
    gsConstantFunction<> gALE(0.,0.,2);
    // boundary conditions: beam
    gsBoundaryConditions<> bcInfoBeam;
    for (index_t d = 0; d < 2; ++d)
        bcInfoBeam.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoALE;
    for (auto it = geoALE.bBegin(); it != geoALE.bEnd(); ++it)
        for (index_t d = 0; d < 2; ++d)
            bcInfoALE.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    gsInfo << "Initialized elasticity system for the beam with " << elAssembler.numDofs() << " dofs.\n";

    index_t numDofsALE = 0;
    for (index_t p = 0; p < geoALE.nPatches(); ++p)
        numDofsALE += geoALE.patch(p).coefsSize()*geoALE.patch(p).coefDim();
    gsInfo << "Initialized elasticity system for the ALE mesh with  " << numDofsALE << " dofs (approx).\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    gsMultiPatch<> displacement,ALE;
    // isogeometric fields (geometry + solution)
    gsField<> displacementField(geoBeam,displacement);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    // paraview collection of time steps
    gsParaviewCollection collectionBeam("flappingBeam_ALE_beam");
    gsParaviewCollection collectionALE("flappingBeam_ALE_mesh");

    gsProgressBar bar;
    gsStopwatch clock;

    // write time and deformation norm
    std::ofstream file;
    file.open("flappingBeam_ALE.txt");
    file << 0. << " " << 0. << std::endl;

    //=============================================//
             // Initial conditions //
    //=============================================//

    // plotting initial condition
    elAssembler.constructSolution(gsMatrix<>::Zero(elAssembler.numDofs(),1),displacement);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_ALE_beam",collectionBeam,0,numPlotPoints);
        plotGeometry(geoALE,"flappingBeam_ALE_mesh",collectionALE,0);
    }

    std::vector<gsMatrix<> > interfaceNow, interfaceNew;
    interfaceNow.push_back(displacement.patch(0).boundary(boundary::north)->coefs());
    interfaceNow.push_back(displacement.patch(0).boundary(boundary::south)->coefs());
    interfaceNow.push_back(displacement.patch(0).boundary(boundary::east)->coefs());
    interfaceNew = interfaceNow;

    elTimeSolver.initialize();

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    clock.restart();
    gsInfo << "Running...\n";
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));
        index_t bad = checkGeometry(geoALE);
        if (bad != -1 && stop)
        {
            gsInfo << "\n Bad patch: " << bad << std::endl;
            break;
        }

        // BEAM
        elTimeSolver.makeTimeStep(timeStep);
        elAssembler.constructSolution(elTimeSolver.displacementVector(),displacement);

        interfaceNew[0] = displacement.patch(0).boundary(boundary::north)->coefs();
        interfaceNew[1] = displacement.patch(0).boundary(boundary::south)->coefs();
        interfaceNew[2] = displacement.patch(0).boundary(boundary::east)->coefs();

        // ALE
        gsElasticityAssembler<real_t> aleAssembler(geoALE,basisALE,bcInfoALE,gALE);
        aleAssembler.options().setReal("PoissonsRatio",poissonsRatioMesh);
        aleAssembler.options().setReal("LocalStiff",stiffDegree);
        aleAssembler.setFixedDofs(3,boundary::south,interfaceNew[0]-interfaceNow[0]);
        aleAssembler.setFixedDofs(4,boundary::north,interfaceNew[1]-interfaceNow[1]);
        aleAssembler.setFixedDofs(5,boundary::west,interfaceNew[2]-interfaceNow[2]);
        aleAssembler.assemble();
        gsSparseSolver<>::PardisoLU solverALE(aleAssembler.matrix());
        gsVector<> solVector = solverALE.solve(aleAssembler.rhs());
        // construct ALE update
        gsMultiPatch<> aleUpdate;
        aleAssembler.constructSolution(solVector,aleAssembler.allFixedDofs(),aleUpdate);
        if (i == 0)
            aleAssembler.constructSolution(solVector,aleAssembler.allFixedDofs(),ALE);
        // apply new deformation to the ALE domain
        for (index_t p = 0; p < geoALE.nPatches(); ++p)
        {
            if (i > 0)
                ALE.patch(p).coefs() += aleUpdate.patch(p).coefs();
            geoALE.patch(p).coefs() += aleUpdate.patch(p).coefs();
        }

        interfaceNow = interfaceNew;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_ALE_beam",collectionBeam,i+1,numPlotPoints);
            plotGeometry(geoALE,"flappingBeam_ALE_mesh",collectionALE,i+1);
        }

        real_t aleNorm = 0.;
        for (index_t p = 0; p < ALE.nPatches(); ++p)
            aleNorm += pow(ALE.patch(p).coefs().norm(),2);

        file << timeStep*(i+1) << " " << sqrt(aleNorm) << std::endl;
    }

    gsInfo << "Complete in " << clock.stop() << "s.\n";
    //=============================================//
                   // Final touches //
    //=============================================//

    if (numPlotPoints > 0)
    {
        gsInfo << "Open \"flappingBeam_mesh_*.pvd\" in Paraview for visualization.\n";
        collectionBeam.save();
        collectionALE.save();
    }

    file.close();

    return 0;
}
