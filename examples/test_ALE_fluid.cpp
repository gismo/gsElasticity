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
#include <gsElasticity/gsALE2.h>
#include <gsElasticity/gsPartitionedFSI2.h>
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
    index_t ALEmethod = ale_method::TEST;
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
    
    gsMultiPatch<> geoALE;
    gsMultiPatch<> geoSquare;
    geoSquare.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,0));

	gsMatrix<> dispBeam = geoSquare.patch(0).coefs();


    for (index_t p = 0; p < 3; ++p)
        geoALE.addPatch(geoFlow.patch(p+3).clone());
    geoALE.computeTopology();

    for (index_t i = 0; i < numUniRef; ++i)
    {
        // basisDisplacement.uniformRefine();
        geoFlow.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisPressure(geoFlow);
    // I use subgrid elements (because degree elevation is not implemented for geometries, so I cant use Taylor-Hood)
    // basisDisplacement.uniformRefine();
    geoALE.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisALE(geoALE);

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
    
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(0,boundary::south,condition_type::dirichlet,0,d);

        bcInfoFlow.addCondition(1,boundary::west,condition_type::dirichlet,0,d);


        bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::west,condition_type::dirichlet,0,d);


        bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);


        bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);

        bcInfoFlow.addCondition(5,boundary::east,condition_type::dirichlet,0,d);        
        bcInfoFlow.addCondition(5,boundary::south,condition_type::dirichlet,0,d);


        bcInfoFlow.addCondition(6,boundary::east,condition_type::dirichlet,0,d);

        bcInfoFlow.addCondition(7,boundary::east,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(7,boundary::north,condition_type::dirichlet,0,d);

    }

    // ALE to flow bdry interface: Navier-Stokes solver contatains a reference to the ALE velocity field,
    // and the FSI module use the ALE displacement to deform the flow geometry
    
    //This might be incorrect, the interface is Structure--ALE--Fluid, so it should be another side of the ALE

    gsBoundaryInterface interfaceALE2Flow;
    interfaceALE2Flow.addInterfaceSide(1,boundary::east,1,boundary::east);
    interfaceALE2Flow.addInterfaceSide(3,boundary::north,3,boundary::north);
    interfaceALE2Flow.addInterfaceSide(4,boundary::south,4,boundary::south);
    interfaceALE2Flow.addInterfaceSide(6,boundary::west,6,boundary::west);



    // interfaceALE2Flow.addPatches(0,3);
    // interfaceALE2Flow.addPatches(1,4);
    // interfaceALE2Flow.addPatches(2,5);

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
    nsTimeSolver.options().setSwitch("ALE",true);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
  
    // mesh deformation module


    // gsALE2<real_t> moduleALE(geoALE,dispBeam,ale_method::method(ALEmethod));
    // moduleALE.options().setReal("LocalStiff",meshStiff);
    // moduleALE.options().setReal("PoissonsRatio",meshPR);
    // moduleALE.options().setSwitch("Check",false);
    // gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";
    // // FSI coupling module
    // gsPartitionedFSI2<real_t> moduleFSI(nsTimeSolver,velFlow, presFlow,moduleALE,dispALE,velALE);

    // moduleFSI.options().setInt("MaxIter",maxCouplingIter);
    // moduleFSI.options().setReal("AbsTol",1e-10);
    // moduleFSI.options().setReal("RelTol",1e-6);
    // moduleFSI.options().setInt("Verbosity",verbosity);



  }