/// This is an example of solving the steady 2D Fluid-Structure Interaction problem
/// using the incompressible Navier-Stokes and nonlinear elasticity solvers
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

void refineBoundaryLayer(gsMultiBasis<> & velocity, gsMultiBasis<> & pressure)
{
    gsMatrix<> boxEast(2,2);
    boxEast << 0.8,1.,0.,0.;
    velocity.refine(0,boxEast);
    pressure.refine(0,boxEast);

    gsMatrix<> boxSouth(2,2);
    boxSouth << 0.,0.,0.,0.2;
    velocity.refine(1,boxSouth);
    velocity.refine(3,boxSouth);
    pressure.refine(1,boxSouth);
    pressure.refine(3,boxSouth);

    gsMatrix<> boxNorth(2,2);
    boxNorth << 0.,0.,0.8,1.;
    velocity.refine(2,boxNorth);
    velocity.refine(4,boxNorth);
    pressure.refine(2,boxNorth);
    pressure.refine(4,boxNorth);

    gsMatrix<> boxWest(2,2);
    boxWest << 0.,0.2,0.,0.;
    velocity.refine(5,boxWest);
    pressure.refine(5,boxWest);
}

int main(int argc, char* argv[]){
    gsInfo << "Testing the steady fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filenameFlow = ELAST_DATA_DIR"/fsi_flow_around_cylinder.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/fsi_beam_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements
    index_t numPlotPoints = 10000;
    real_t youngsModulus = 0.5e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    bool supg = false;
    index_t iter = 5;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the time-dependent Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("b","blayer","Number of additional boundary layer refinements",numBLRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("y","young","Young's modulus of the beam materail",youngsModulus);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
    cmd.addInt("i","iter","Number of coupling iterations",iter);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

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

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoMesh;
    for (auto it = geoFlow.bBegin(); it != geoFlow.bEnd(); ++it)
        for (index_t d = 0; d < 2; ++d)
            bcInfoMesh.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    // creating bases
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisPressure(geoFlow);
    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisDisplacement.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);
    // additional velocity refinement for stable mixed FEM;
    // displacement is also refined to keep bases match at the fsi interface
    if (subgrid)
    {
        basisVelocity.uniformRefine();
        basisDisplacement.uniformRefine();
    }
    else
    {
        basisVelocity.degreeElevate();
        basisDisplacement.degreeElevate();
    }
    // navier stokes assembler
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,g);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: flow mesh
    gsElasticityAssembler<real_t> meshAssembler(geoFlow,basisVelocity,bcInfoMesh,g);
    elAssembler.options().setReal("PoissonsRatio",0.4);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsInfo << "Initialized elasticity system for mesh deformation with " << meshAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // containers for solution as IGA functions
    gsMultiPatch<> velocity, pressure, displacementBeam, displacementMesh;
    // isogeometric fields (geometry + solution)
    gsField<> velocityField(geoFlow,velocity);
    gsField<> pressureField(geoFlow,pressure);
    gsField<> displacementFieldBeam(geoBeam,displacementBeam);
    gsField<> displacementFieldMesh(geoFlow,displacementMesh);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    fieldsFlow["Displacement"] = &displacementFieldMesh;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementFieldBeam;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("fsi_steady_flow");
    gsParaviewCollection collectionBeam("fsi_steady_beam");
    // plotting initial condition
    nsAssembler.constructSolution(gsMatrix<>::Zero(nsAssembler.numDofs(),1),velocity,pressure);
    elAssembler.constructSolution(gsMatrix<>::Zero(elAssembler.numDofs(),1),displacementBeam);
    meshAssembler.constructSolution(gsMatrix<>::Zero(meshAssembler.numDofs(),1),displacementMesh);
    gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"fsi_steady_flow",collectionFlow,0,numPlotPoints);
    gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"fsi_steady_beam",collectionBeam,0,numPlotPoints);

    //=============================================//
             // Coupled simulation //
    //=============================================//

    gsStopwatch clock;
    clock.restart();
    for (index_t i = 0; i < iter; ++i)
    {
        // 1. apply flow-mesh deformation (how?)
        // 2. solve flow (done)
        // 3. compute drag&lift (done)
        // 4. apply drag&lift (how?)
        // 5. solve beam (done)
        // 6. compute flow mesh deformation (done)
    }
    gsInfo << "Solverd in "<< clock.stop() <<"s.\n";

    //=============================================//
             // Final touches //
    //=============================================//

    gsInfo << "Plotting the output to the Paraview files \"fsi_steady flow.pvd\" and \"fsi_steady beam.pvd\"...\n";
    collectionFlow.save();
    collectionBeam.save();

    return 0;
}
