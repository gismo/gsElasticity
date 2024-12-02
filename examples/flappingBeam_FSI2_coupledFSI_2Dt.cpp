/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsALE.h>
#include <gsElasticity/gsPartitionedFSI.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

void writeLog(std::ofstream & ofs, const gsNsAssembler<real_t> & assemblerFlow,
              const gsMultiPatch<> & velocity, const gsMultiPatch<> & pressure,
              const gsMultiPatch<> & displacementBeam, const gsMultiPatch<> & geoALE, const gsMultiPatch<> & dispALE,
              real_t simTime, real_t aleTime, real_t flowTime, real_t beamTime,
              index_t couplingIter, index_t flowIter, index_t beamIter,
              real_t omega, real_t resAbs, real_t resRel)
{
    // computing force acting on the surface of the submerged structure
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = assemblerFlow.computeForce(velocity,pressure,bdrySides);

    // compute the pressure difference between the front and the end points of the structure
    gsMatrix<> presFront(2,1);
    presFront << 1.,0.5;
    presFront = pressure.patch(0).eval(presFront);
    gsMatrix<> presBack(2,1);
    presBack << 0.,0.5;
    presBack = pressure.patch(5).eval(presBack);

    // compute displacement of the beam point A
    gsMatrix<> dispA(2,1);
    dispA << 1.,0.5;
    dispA = displacementBeam.patch(0).eval(dispA);

    // print: 1-simTime 2-drag 3-lift 4-pressureDiff 5-dispAx 6-dispAy 7-aleNorm
    //        8-aleTime 9-flowTime 10-beamTime 11-couplingIter 12-flowIter 13-beamIter
    //        14-omega 15-resAbs 16-resRel
    ofs << simTime << " " << force.at(0) << " " << force.at(1) << " " << presFront.at(0)-presBack.at(0) << " "
        << dispA.at(0) << " " << dispA.at(1) << " " << normL2(geoALE,dispALE) << " "
        << aleTime << " " << flowTime << " " << beamTime << " "
        << couplingIter << " " << flowIter << " " << beamIter << " "
        << omega << " " << resAbs << " " << resRel << std::endl;
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the two-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // beam parameters
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t densitySolid = 1.0e4;
    real_t beamLoad = 0.0;
    // flow parameters
    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flow.xml";
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    // ALE parameters
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.5;
    index_t ALEmethod = ale_method::LE;
    bool oneWay = false;
    // space discretization
    index_t numUniRef = 4;
    // time integration
    real_t timeStep = 0.0025;
    real_t timeSpan = 15.;
    real_t thetaFluid = 0.5;
    real_t thetaSolid = 1.;
    index_t maxCouplingIter = 10;
    bool imexOrNewton = true;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 0.;
    index_t verbosity = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("l","load","Gravitation acceleration acting on the beam",beamLoad);
    cmd.addReal("x","chi","Local stiffening degree for ALE",meshStiff);
    cmd.addSwitch("o","oneway","Run as a oneway coupled simulation: beam-to-flow",oneWay);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addInt("i","iter","Number of coupling iterations",maxCouplingIter);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addInt("v","verbosity","Amount of info printed to the prompt: 0 - none, 1 - crucial, 2 - all",verbosity);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (oneWay)
    {
        maxCouplingIter = 1;
        densitySolid = 1.0e3;
        beamLoad = 2*densitySolid;
        thetaSolid = 0.5;
    }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);
    // I deform only three flow patches adjacent to the FSI interface
    gsMultiPatch<> geoALE;
    for (index_t p = 0; p < 3; ++p)
        geoALE.addPatch(geoFlow.patch(p+3).clone());
    geoALE.computeTopology();


    // gsWriteParaview(geoALE,"geoALE",1000);
    // return 0;

    // creating bases
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        geoFlow.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisPressure(geoFlow);
    // I use subgrid elements (because degree elevation is not implemented for geometries, so I cant use Taylor-Hood)
    basisDisplacement.uniformRefine();
    geoALE.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> gZero(0.,0.,2);
    gsConstantFunction<> g(0.,beamLoad,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velFlow, presFlow, dispBeam, dispALE, velALE;
    // boundary conditions: flow
    gsBoundaryConditions<> bcInfoFlow;
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0);
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfoFlow.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(5,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::north,condition_type::dirichlet,0,d);
    }
    // boundary conditions: beam
    gsBoundaryConditions<> bcInfoBeam;
    for (index_t d = 0; d < 2; ++d)
        bcInfoBeam.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
    // flow to beam interface: these Neumann boundary condtions contain references to flow and ALE solutions;
    // by updating them, we update the boundary load as well
    gsFsiLoad<real_t> fSouth(geoALE,dispALE,1,boundary::north,
                             velFlow,presFlow,4,viscosity,densityFluid);
    gsFsiLoad<real_t> fEast(geoALE,dispALE,2,boundary::west,
                            velFlow,presFlow,5,viscosity,densityFluid);
    gsFsiLoad<real_t> fNorth(geoALE,dispALE,0,boundary::south,
                             velFlow,presFlow,3,viscosity,densityFluid);
    if (!oneWay)
    {
        bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&fSouth);
        bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&fEast);
        bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&fNorth);
    }

    // beam to ALE interface: ALE module contains a reference to the beam displacement field;
    // by updating the displacement field, we update the displacement of the FSI interface in ALE computations
    gsBoundaryInterface interfaceBeam2ALE;
    interfaceBeam2ALE.addInterfaceSide(0,boundary::north,0,boundary::south);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::south,1,boundary::north);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::east,2,boundary::west);

    // ALE to flow bdry interface: Navier-Stokes solver contatains a reference to the ALE velocity field,
    // and the FSI module use the ALE displacement to deform the flow geometry
    gsBoundaryInterface interfaceALE2Flow;
    interfaceALE2Flow.addInterfaceSide(0,boundary::south,3,boundary::south);
    interfaceALE2Flow.addInterfaceSide(1,boundary::north,4,boundary::north);
    interfaceALE2Flow.addInterfaceSide(2,boundary::west,5,boundary::west);
    interfaceALE2Flow.addPatches(0,3);
    interfaceALE2Flow.addPatches(1,4);
    interfaceALE2Flow.addPatches(2,5);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,gZero);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,gZero);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler,&velALE,&interfaceALE2Flow);
    nsTimeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    nsTimeSolver.options().setReal("Theta",thetaFluid);
    nsTimeSolver.options().setSwitch("ALE",true);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gZero);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    elTimeSolver.options().setReal("Beta",thetaSolid/2);
    elTimeSolver.options().setReal("Gamma",thetaSolid);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // mesh deformation module
    gsALE<real_t> moduleALE(geoALE,dispBeam,interfaceBeam2ALE,ale_method::method(ALEmethod));
    moduleALE.options().setReal("LocalStiff",meshStiff);
    moduleALE.options().setReal("PoissonsRatio",meshPR);
    moduleALE.options().setSwitch("Check",false);
    gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";
    // FSI coupling module
    gsPartitionedFSI<real_t> moduleFSI(nsTimeSolver,velFlow, presFlow,
                                       elTimeSolver,dispBeam,
                                       moduleALE,dispALE,velALE);
    moduleFSI.options().setInt("MaxIter",maxCouplingIter);
    moduleFSI.options().setReal("AbsTol",1e-10);
    moduleFSI.options().setReal("RelTol",1e-6);
    moduleFSI.options().setInt("Verbosity",verbosity);

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // beam stress field
    gsPiecewiseFunction<> stresses;
    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velFlow);
    gsField<> pressureField(nsAssembler.patches(),presFlow);
    gsField<> displacementField(geoBeam,dispBeam);
    gsField<> stressField(geoBeam,stresses,true);
    gsField<> aleField(geoALE,dispALE);

    // creating containers to plot several fields corresponding to the same geometry to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    fieldsBeam["von Mises"] = &stressField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSI2_flow");
    gsParaviewCollection collectionBeam("flappingBeam_FSI2_beam");
    gsParaviewCollection collectionALE("flappingBeam_FSI2_ALE");

    std::ofstream logFile;
    logFile.open("flappingBeam_FSI2.txt");
    logFile << "# simTime drag lift presDiff dispAx dispAy aleNorm aleTime flowTime"
            << " beamTime couplingIter flowIter beamIter omega resAbs resRel\n";

    gsProgressBar bar;
    gsStopwatch totalClock, iterClock;

    //=============================================//
                   // Initial condtions //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.getFixedDofs(0,boundary::west,inflowDDoFs);
    nsAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());

    elTimeSolver.setDisplacementVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));
    elTimeSolver.setVelocityVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));

    // plotting initial condition
    nsTimeSolver.constructSolution(velFlow,presFlow);
    elTimeSolver.constructSolution(dispBeam);
    elAssembler.constructCauchyStresses(dispBeam,stresses,stress_components::von_mises);
    moduleALE.constructSolution(dispALE);

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,0,numPlotPoints);
        plotDeformation(geoALE,dispALE,"flappingBeam_FSI2_ALE",collectionALE,0);
    }
    writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,0.,0.,0.,0.,0,0,0,1.,0.,0.);

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    real_t simTime = 0.;
    real_t numTimeStep = 0;
    real_t timeALE = 0.;
    real_t timeFlow = 0.;
    real_t timeBeam = 0.;

    totalClock.restart();

    gsInfo << "Running the simulation...\n";
    while (simTime < timeSpan)
    {
        bar.display(simTime/timeSpan);
        // change time step for the initial warm-up phase
        real_t tStep = (warmUp && simTime < 2.) ? 0.1 : timeStep;
        // smoothly change the inflow boundary condition
        if (simTime < 2.)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(EIGEN_PI*(simTime+tStep)/2.))/2);
        if (simTime > 7.)
            moduleALE.options().setSwitch("Check",true);

        if (!moduleFSI.makeTimeStep(tStep))
        {
            gsInfo << "Invalid ALE mapping. Terminated.\n";
            break;
        }

        // Iteration end
        simTime += tStep;
        timeALE += moduleFSI.timeALE();
        timeBeam += moduleFSI.timeEL();
        timeFlow += moduleFSI.timeNS();
        numTimeStep++;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,numTimeStep,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_FSI2_beam",collectionBeam,numTimeStep,numPlotPoints);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_FSI2_ALE",collectionALE,numTimeStep,numPlotPoints);
            plotDeformation(geoALE,dispALE,"flappingBeam_FSI2_ALE",collectionALE,numTimeStep);
        }
        writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,
                 simTime,timeALE,timeFlow,timeBeam, moduleFSI.numberIterations(),
                 nsTimeSolver.numberIterations(),elTimeSolver.numberIterations(),
                 moduleFSI.aitkenOmega(),moduleFSI.residualNormAbs(),moduleFSI.residualNormRel());
    }

    //=============================================//
                   // Final touches //
    //=============================================//

    gsInfo << "Complete in: " << secToHMS(totalClock.stop())
           << ", ALE time: " << secToHMS(timeALE)
           << ", flow time: " << secToHMS(timeFlow)
           << ", beam time: " << secToHMS(timeBeam) << std::endl;

    if (numPlotPoints > 0)
    {
        collectionFlow.save();
        collectionBeam.save();
        collectionALE.save();
        gsInfo << "Open \"flappingBeam_FSI2_*.pvd\" in Paraview for visualization.\n";
    }
    logFile.close();
    gsInfo << "Log file created in \"flappingBeam_FSI2.txt\".\n";

    return 0;
}
