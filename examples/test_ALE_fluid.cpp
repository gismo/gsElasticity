/**
 *
 *
 *
 *
 *
*/

#include <gismo.h>
//#include <gsElasticity/gsElasticityAssembler.h>
//#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsALE.h>
// #include <gsElasticity/gsPartitionedFSI2.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;




int main(int argc, char* argv[])
{
    gsInfo << "Testing the two-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // beam parameters
    // std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    // real_t youngsModulus = 1.4e6;
    // real_t poissonsRatio = 0.4;
    // real_t densitySolid = 1.0e4;
    // flow parameters
    std::string filenameFlow = ELAST_DATA_DIR"/grid.xml";
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    // ALE parameters
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.5;
    index_t ALEmethod = ale_method::TINE;
    bool oneWay = false;
    // space discretization
    index_t numUniRef = 3;
    // time integration
    real_t timeStep = 0.01;
    real_t timeSpan = 15.;
    real_t thetaFluid = 0.5;
    index_t maxCouplingIter = 10;
    bool imexOrNewton = false;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 0.;
    index_t verbosity = 0;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("x","chi","Local stiffening degree for ALE",meshStiff);
    cmd.addSwitch("o","oneway","Run as a oneway coupled simulation: beam-to-flow",oneWay);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE, 6 - TEST",ALEmethod);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addInt("i","iter","Number of coupling iterations",maxCouplingIter);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addInt("v","verbosity","Amount of info printed to the prompt: 0 - none, 1 - crucial, 2 - all",verbosity);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//
    // Scanning geometry
	gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);    

    // creating bases
    // gsMultiBasis<> basisDisplacement(geoBeam);
    
    gsMultiPatch<> geoALE = geoFlow;
    gsMultiPatch<> geoSquare;
    geoSquare.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,0));

	// gsMatrix<> dispBeam = geoSquare.patch(0).coefs();

    // for (index_t p = 0; p < 3; ++p)
    //     geoALE.addPatch(geoFlow.patch(p+3).clone());
    // geoALE.computeTopology();

    for (index_t i = 0; i < numUniRef; ++i)
    {
        // basisDisplacement.uniformRefine();
        geoFlow.uniformRefine();
        geoSquare.uniformRefine();
        geoALE.uniformRefine();
    }

    gsMultiPatch<> dispBeam = geoSquare;

    gsMultiBasis<> basisPressure(geoFlow);
    // I use subgrid elements (because degree elevation is not implemented for geometries, so I cant use Taylor-Hood)
    // basisDisplacement.uniformRefine();

    gsWarn<<"Fix the thing below this line\n";
    // geoALE.uniformRefine();
    // geoSquare.uniformRefine();
    // geoFlow.uniformRefine();
        geoALE.degreeElevate();
        geoSquare.degreeElevate();
        geoFlow.degreeElevate();

    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> gZero(0.,0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    // gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velFlow, presFlow, dispALE, velALE;
    // boundary conditions: flow
    gsBoundaryConditions<> bcInfoFlow;
    // bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
    // bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    

    /*
     * ______________________
     * |      |      |      |
     * |  2   |  4   |  7   |
     * |______|______|______|
     * |      |      |      |
     * |  1   |      |  6   |
     * |______|______|______|
     * |      |      |      |
     * |  0   |  3   |  5   |
     * |______|______|______|
     *
     */
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        // bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
        // bcInfoFlow.addCondition(0,boundary::south,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(1,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::east,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        // bcInfoFlow.addCondition(2,boundary::west,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::north,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::south,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(5,boundary::east,condition_type::dirichlet,0,d);
        // bcInfoFlow.addCondition(5,boundary::south,condition_type::dirichlet,0,d);

        bcInfoFlow.addCondition(6,boundary::west,condition_type::dirichlet,0,d);
        // bcInfoFlow.addCondition(6,boundary::east,condition_type::dirichlet,0,d);

        // bcInfoFlow.addCondition(7,boundary::east,condition_type::dirichlet,0,d);
        // bcInfoFlow.addCondition(7,boundary::north,condition_type::dirichlet,0,d);
    }

    // ALE to flow bdry interface: Navier-Stokes solver contatains a reference to the ALE velocity field,
    // and the FSI module use the ALE displacement to deform the flow geometry
    
    //This might be incorrect, the interface is Structure--ALE--Fluid, so it should be another side of the ALE

    gsBoundaryInterface interfaceALE2Flow;
    interfaceALE2Flow.addInterfaceSide(1,boundary::east,1,boundary::east);
    interfaceALE2Flow.addInterfaceSide(3,boundary::north,3,boundary::north);
    interfaceALE2Flow.addInterfaceSide(4,boundary::south,4,boundary::south);
    interfaceALE2Flow.addInterfaceSide(6,boundary::west,6,boundary::west);


    gsBoundaryInterface interfaceSolid2ALE;
    interfaceSolid2ALE.addInterfaceSide(0,boundary::west,1,boundary::east);
    interfaceSolid2ALE.addInterfaceSide(0,boundary::east,6,boundary::west);
    interfaceSolid2ALE.addInterfaceSide(0,boundary::south,3,boundary::north);
    interfaceSolid2ALE.addInterfaceSide(0,boundary::north,4,boundary::south);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

     // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,gZero);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,gZero);

    nsMassAssembler.options().setReal("Density",densityFluid);
    // velALE: velocity update from the ALE solver
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler,&velALE,&interfaceALE2Flow);

    nsTimeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    nsTimeSolver.options().setReal("Theta",thetaFluid);
    nsTimeSolver.options().setInt("Verbosity",verbosity);
    gsDebugVar(nsTimeSolver.options().getInt("Verbosity"));
    nsTimeSolver.options().setSwitch("ALE",true);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
  
    // mesh deformation module


    gsALE<real_t> moduleALE(geoALE,dispBeam,interfaceSolid2ALE,ale_method::method(ALEmethod));
    moduleALE.options().setReal("LocalStiff",meshStiff);
    moduleALE.options().setReal("PoissonsRatio",meshPR);
    moduleALE.options().setSwitch("Check",false);
    // gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";
    // // FSI coupling module
    // gsPartitionedFSI2<real_t> moduleFSI(nsTimeSolver,velFlow, presFlow,moduleALE,dispALE,velALE);

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velFlow);
    gsField<> pressureField(nsAssembler.patches(),presFlow);
    gsField<> aleField(geoALE,dispALE);

    // creating containers to plot several fields corresponding to the same geometry to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;

    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSI2_flow");
    gsParaviewCollection collectionALE("flappingBeam_FSI2_ALE");


    //=============================================//
                   // Initial condtions //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());

    //=============================================//
                   // Uncoupled simulation //
    //=============================================//

    real_t time = 0;
    real_t angle = 0;
    real_t maxAngle = 20;
    gsMultiPatch<> geoFlowDisp, geoFlowVelo;
    std::string fn = "output";
    gsParaviewCollection collection(fn);

    index_t kmax = math::ceil(timeSpan/timeStep)+1;
    gsMultiPatch<> geoSquareNew = geoSquare;
    gsMultiPatch<> geoSquareOld = geoSquare;
    dispBeam = geoSquare;
    for (index_t k = 0; k!=kmax; k++)
    {
        gsInfo << "Time: " << time << "\n";

        // Rotate geometry
        angle = maxAngle * math::sin(2*M_PI/timeSpan*time);
        // angle = 0.1 * math::sin(2*M_PI/timeSpan*time);
        // gsDebugVar(angle);

        // Displacements of the geometry
        geoSquareNew = geoSquare;
        // gsNurbsCreator<>::shift2D(geoSquareNew.patch(0),angle,0);
        gsNurbsCreator<>::rotate2D(geoSquareNew.patch(0),angle,0.5,0.5);
        dispBeam.patch(0).coefs() = geoSquareNew.patch(0).coefs() - geoSquareOld.patch(0).coefs();

        // Store the old fluid mesh in the velocity mesh
        moduleALE.constructSolution(geoFlowVelo);

        // Compute the fluid mesh displacements
        moduleALE.updateMesh();

        // Update the fluid mesh displacements
        moduleALE.constructSolution(geoFlowDisp);

        // Update the fluid mesh velocities (v = (u_t - u_t-1) / dt)
        for (size_t p = 0; p < geoFlowVelo.nPatches(); ++p)
            geoFlowVelo.patch(p).coefs() = (geoFlowDisp.patch(p).coefs() - geoFlowVelo.patch(p).coefs()) / timeStep;

        // Update the fluid mesh
        for (size_t p = 0; p!=geoFlow.nPatches(); ++p)
            geoFlow.patch(p).coefs() += geoFlowDisp.patch(p).coefs();

        // Update geoALE (FOR TESTING)
        geoALE = geoFlow;
        dispALE = geoFlowDisp;

        ///////// We might need something like "recover state"


        // Update fluid mesh
        // Update the fluid mesh
        for (size_t p = 0; p!=geoFlow.nPatches(); ++p)
        {
            nsTimeSolver.assembler().patches().patch(p).coefs() = geoFlow.patch(p).coefs();
            nsTimeSolver.mAssembler().patches().patch(p).coefs() = geoFlow.patch(p).coefs();
        }


        // set velocity boundary condition on the FSI interface; velocity comes from the ALE velocity;
        // FSI inteface info is contained in the Navier-Stokes solver
        for (size_t p = 0; p < nsTimeSolver.aleInterface().sidesA.size(); ++p)
        {
            index_t pFlow = nsTimeSolver.aleInterface().sidesB[p].patch;
            boxSide sFlow = nsTimeSolver.aleInterface().sidesB[p].side();
            index_t pALE = nsTimeSolver.aleInterface().sidesA[p].patch;
            boxSide sALE = nsTimeSolver.aleInterface().sidesA[p].side();
            nsTimeSolver.assembler().setFixedDofs(pFlow,sFlow,geoFlowVelo.patch(pALE).boundary(sALE)->coefs());
        }


        // Actual solve
        nsTimeSolver.makeTimeStep(timeStep);
        nsTimeSolver.constructSolution(velFlow,presFlow);

        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flappingBeam_FSI2_flow",collectionFlow,k,1000);









        const std::string fileName = fn + util::to_string(k);
        gsWriteParaview(geoFlow,fileName,1000,true);
        for (size_t p = 0; p!=geoFlow.nPatches(); ++p)
        {
            const std::string fileName_nopath = gsFileManager::getFilename(fileName);
            gsWriteParaview(geoFlow.patch(p),geoFlow.patch(p).support(),fileName,1000,true);
            collection.addPart(fileName_nopath + "_" + std::to_string(p) + ".vts",k,"geometry",p);
            collection.addPart(fileName_nopath + "_" + std::to_string(p) + "_mesh" + ".vtp",k,"mesh",p);
        }

        time += timeStep;

        geoSquareOld = geoSquareNew;
    }

    collection.save();
    collectionFlow.save();


    for (size_t p = 0; p!=geoFlow.nPatches(); ++p)
        geoFlow.patch(p).coefs() += geoFlowDisp.patch(p).coefs();



    // moduleFSI.options().setInt("MaxIter",maxCouplingIter);
    // moduleFSI.options().setReal("AbsTol",1e-10);
    // moduleFSI.options().setReal("RelTol",1e-6);
    // moduleFSI.options().setInt("Verbosity",verbosity);



  }