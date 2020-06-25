/// This is the benchmark for computing ALE mapping for an FSI problem based on this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// There is no flow in the benchmark. The flow domain deformation is driven by a freely oscillating beam.
/// The ALE mapping is computed using the nonlinear elasticity method with Jacobian-based local stiffening.
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsALE.h>

using namespace gismo;

void writeLog(std::ofstream & ofs,
              const gsMultiPatch<> & displacementBeam, const gsMultiPatch<> & geoALE, const gsMultiPatch<> & dispALE,
              real_t simTime, real_t aleTime, real_t beamTime, index_t beamIter)
{
    // compute displacement of the beam point A
    gsMatrix<> dispA(2,1);
    dispA << 1.,0.5;
    dispA = displacementBeam.patch(0).eval(dispA);

    // print: simTime dispAx dispAy aleNorm aleTime beamTime beamIter
    ofs << simTime << " " << dispA.at(0) << " " << dispA.at(1) << " " << normL2(geoALE,dispALE)
        << " " << aleTime << " " << beamTime << " " << beamIter << std::endl;
}


int main(int argc, char* argv[])
{
    gsInfo << "Testing the ALE mapping construction in 2D.\n";

    std::string filenameALE = ELAST_DATA_DIR"/flappingBeam_flow.xml";
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
    index_t ALEmethod = ale_method::TINE;
    bool check = true;
    index_t numIter = 1;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the steady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Pois ratio for mesh",poissonsRatioMesh);
    cmd.addReal("l","load","Gravitational loading acting on the beam",loading);
    cmd.addReal("x","xjac","Stiffening degree for the Jacobian-based local stiffening",stiffDegree);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
    cmd.addSwitch("c","check","Check bijectivity of the ALE displacement field",check);
    cmd.addInt("i","iter","Number of iterations for nonlinear methods",numIter);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);
    gsMultiPatch<> geoFlow; // this is a full flow geometry
    gsReadFile<>(filenameALE, geoFlow);
    // this is a copy of the flow domain with only the first 6 patches.
    gsMultiPatch<> geoALE;
    for (index_t p = 0; p < 3; ++p)
        geoALE.addPatch(geoFlow.patch(p+3).clone());
    geoALE.computeTopology();

    // creating bases
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsMultiPatch<> displacement,ALE;
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

    gsBoundaryInterface interfaceBeam2ALE;
    interfaceBeam2ALE.addInterfaceSide(0,boundary::north,0,boundary::south);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::south,1,boundary::north);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::east,2,boundary::west);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gBeam);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";

    // ALE module with nonlinear elasticity
    gsALE<real_t> moduleALE(geoALE,displacement,interfaceBeam2ALE,ale_method::method(ALEmethod));
    moduleALE.options().setReal("LocalStiff",stiffDegree);
    moduleALE.options().setReal("PoissonsRatio",poissonsRatioMesh);
    moduleALE.options().setSwitch("Check",check);
    moduleALE.options().setInt("NumIter",numIter);
    gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> displacementField(geoBeam,displacement);
    gsField<> aleField(geoALE,ALE);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    std::map<std::string,const gsField<> *> fieldsAle;
    fieldsAle["ALE displacement"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionBeam("flappingBeam_ALE_beam");
    gsParaviewCollection collectionMesh("flappingBeam_ALE_mesh");
    gsParaviewCollection collectionALE("flappingBeam_ALE_ale");

    gsProgressBar bar;
    gsStopwatch iterClock, totalClock;

    // write time and deformation norm
    std::ofstream logFile;
    logFile.open("flappingBeam_ALE.txt");
    logFile << "# simTime dispAx dispAy aleNorm aleTime beamTime beamIter \n";

    //=============================================//
             // Initial conditions //
    //=============================================//

    // set initial conditions for the beam
    elTimeSolver.setDisplacementVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));
    elTimeSolver.setVelocityVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));
    // constructing initial fields
    elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
    moduleALE.constructSolution(ALE);
    // plotting initial condition
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_ALE_beam",collectionBeam,0,numPlotPoints);
        plotDeformation(geoALE,ALE,"flappingBeam_ALE_mesh",collectionMesh,0);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsAle,"flappingBeam_ALE_ale",collectionALE,0,numPlotPoints);
    }
    writeLog(logFile,displacement,geoALE,ALE,0.,0.,0.,0);

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    real_t timeALE = 0.;
    real_t timeBeam = 0.;

    gsInfo << "Running the simulation...\n";
    totalClock.restart();
    for (index_t i = 0; i < (index_t)(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,(index_t)(timeSpan/timeStep));

        // BEAM
        iterClock.restart();
        elTimeSolver.makeTimeStep(timeStep);
        elAssembler.constructSolution(elTimeSolver.displacementVector(),elTimeSolver.allFixedDofs(),displacement);
        timeBeam += iterClock.stop();

        // ALE
        iterClock.restart();
        index_t badPatch = moduleALE.updateMesh();
        moduleALE.constructSolution(ALE);
        timeALE += iterClock.stop();

        // output
        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_ALE_beam",collectionBeam,i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_ALE_mesh",collectionMesh,i+1);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsAle,"flappingBeam_ALE_ale",collectionALE,i+1,numPlotPoints);
        }
        writeLog(logFile,displacement,geoALE,ALE,timeStep*(i+1),timeALE,timeBeam,elTimeSolver.numberIterations());

        // test if any patch is not bijective
        if (badPatch != -1)
        {
            gsInfo << "\n Bad patch: " << badPatch << std::endl;
            break;
        }
    }

    //=============================================//
                   // Final touches //
    //=============================================//

    gsInfo << "Complete in: " << secToHMS(totalClock.stop())
           << ", ALE time: " << secToHMS(timeALE)
           << ", beam time: " << secToHMS(timeBeam) << std::endl;

    if (numPlotPoints > 0)
    {
        collectionBeam.save();
        collectionMesh.save();
        collectionALE.save();
        gsInfo << "Open \"flappingBeam_ALE_*.pvd\" in Paraview for visualization.\n";
    }
    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_ALE.txt\".\n";

    return 0;
}
