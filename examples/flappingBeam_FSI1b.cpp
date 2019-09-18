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

void formResidual(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow,
                  gsMatrix<> & residual)
{
    index_t numBdry = interfaceOld.size();
    short_t dim = interfaceOld[0].cols();
    index_t numDoFs = 0;
    for (index_t i = 0; i < numBdry; ++i)
        numDoFs += dim*interfaceOld[i].rows();
    residual.setZero(numDoFs,1);

    index_t index = 0;
    for (index_t i = 0; i < numBdry; ++i)
    {
        index_t size = interfaceOld[i].rows();
        for (index_t d = 0; d < dim; ++d)
        {
            residual.middleRows(index,size) = interfaceNow[i].col(d) - interfaceOld[i].col(d);
            index += size;
        }
    }
}

void aitkenRelaxation(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow,
                      std::vector<gsMatrix<> > & interfaceNew, real_t & omega)
{
    gsMatrix<> residualOld,residualNew;
    formResidual(interfaceOld,interfaceNow,residualOld);
    formResidual(interfaceNow,interfaceNew,residualNew);

    omega = -1.*omega*(residualOld.transpose()*(residualNew-residualOld))(0,0) /
            pow((residualNew-residualOld).norm(),2);

    for (index_t i = 0; i < index_t(interfaceOld.size()); ++i)
    {
        interfaceOld[i] = interfaceNow[i];
        interfaceNow[i] += omega*(interfaceNew[i]-interfaceNow[i]);
    }
}

real_t computeResidual(std::vector<gsMatrix<> > & interfaceOld, std::vector<gsMatrix<> > & interfaceNow)
{
    gsMatrix<> residual;
    formResidual(interfaceOld,interfaceNow,residual);
    return residual.norm();
}

int main(int argc, char* argv[])
{
    gsInfo << "Testing the unsteady fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filenameFlow = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameFlowPart = ELAST_DATA_DIR"/flappingBeam_flowPart.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRefFlow = 3; // number of h-refinements for the fluid
    index_t numKRefFlow = 0; // number of k-refinements for the fluid
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numUniRefBeam = 3; // number of h-refinements for the beam and the ALE mapping
    index_t numKRefBeam = 0; // number of k-refinements for the beam and the ALE mapping
    index_t numPlotPoints = 10000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    bool supg = false;
    real_t densityFluid = 1000.;
    real_t densitySolid = 1000.;
    real_t absTol = 1e-6;
    real_t relTol = 1e-6;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the unsteady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRefFlow);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("b","beamrefine","Number of uniform refinement applications for the beam and ALE",numUniRefBeam);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("y","young","Young's modulus of the beam materail",youngsModulus);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addSwitch("g","supg","Use SUPG stabilization",supg);
    cmd.addReal("d","density","Density of the solid",densitySolid);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoPart; // this is a part of the flow geometry; we deform only this part to save memory and time
    gsReadFile<>(filenameFlowPart, geoPart);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velocity, pressure, displacement, ALE;
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
    gsFsiLoad<real_t> fSouth(geoPart,ALE,1,boundary::north,
                             velocity,pressure,4,viscosity,densityFluid);
    gsFsiLoad<real_t> fEast(geoPart,ALE,2,boundary::west,
                            velocity,pressure,5,viscosity,densityFluid);
    gsFsiLoad<real_t> fNorth(geoPart,ALE,0,boundary::south,
                             velocity,pressure,3,viscosity,densityFluid);
    bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&fSouth);
    bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&fEast);
    bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&fNorth);
    // boundary conditions: flow mesh, set zero dirichlet on the entire boundary
    gsBoundaryConditions<> bcInfoALE;
    for (auto it = geoPart.bBegin(); it != geoPart.bEnd(); ++it)
        for (index_t d = 0; d < 2; ++d)
            bcInfoALE.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    // creating bases
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisPressure(geoFlow);
    for (index_t i = 0; i < numKRefFlow; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRefFlow; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    for (index_t i = 0; i < numBLRef; ++i)
        refineBoundaryLayer(basisVelocity,basisPressure);

    // additional velocity refinement for stable mixed FEM;
    // displacement is also refined to keep bases match at the fsi interface
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    gsMultiBasis<> basisDisplacement(geoBeam);
    for (index_t i = 0; i < numKRefBeam; ++i)
    {
        basisDisplacement.degreeElevate();
        geoPart.degreeElevate();
        geoFlow.degreeElevate();
        basisDisplacement.uniformRefine();
        geoPart.uniformRefine();
        geoFlow.uniformRefine();
    }
    for (index_t i = 0; i < numUniRefFlow; ++i)
    {
        basisDisplacement.uniformRefine();
        geoPart.uniformRefine();
        geoFlow.uniformRefine();
    }
    gsMultiBasis<> basisALE(geoPart);

    // navier stokes assembler in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,g);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    gsMatrix<> solutionFlow = gsMatrix<>::Zero(nsAssembler.numDofs(),1);
    // elasticity assembler: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    gsMatrix<> solutionBeam = gsMatrix<>::Zero(elAssembler.numDofs(),1);
    // elasticity assembler: flow mesh
    gsElasticityAssembler<real_t> aleAssembler(geoPart,basisALE,bcInfoALE,g);
    aleAssembler.options().setReal("PoissonsRatio",0.4);
    aleAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsInfo << "Initialized elasticity system for ALE with " << aleAssembler.numDofs() << " dofs.\n";
    gsMatrix<> solutionALE = gsMatrix<>::Zero(aleAssembler.numDofs(),1);

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velocity);
    gsField<> pressureField(nsAssembler.patches(),pressure);
    gsField<> displacementField(geoBeam,displacement);
    gsField<> aleField(geoPart,ALE);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    std::map<std::string,const gsField<> *> fieldsPart;
    fieldsPart["ALE"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("fsi_steady_flow");
    gsParaviewCollection collectionBeam("fsi_steady_beam");
    gsParaviewCollection collectionFlowPart("fsi_steady_flow_part");
    // plotting initial condition
    nsAssembler.constructSolution(solutionFlow,velocity,pressure);
    elAssembler.constructSolution(solutionBeam,displacement);
    aleAssembler.constructSolution(solutionALE,ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"fsi_steady_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"fsi_steady_beam",collectionBeam,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsPart,"fsi_steady_flow_part",collectionFlowPart,0,numPlotPoints);
    }
    //=============================================//
             // Coupled simulation //
    //=============================================//


    // containers for interface displacement DoFs
    std::vector<gsMatrix<> > interfaceOld, interfaceNow, interfaceNew;
    // push zeros to the old interface
    interfaceOld.push_back(displacement.patch(0).boundary(boundary::north)->coefs());
    interfaceOld.push_back(displacement.patch(0).boundary(boundary::south)->coefs());
    interfaceOld.push_back(displacement.patch(0).boundary(boundary::east)->coefs());
    // Aitken relaxation factor
    real_t omega;
    // convergence
    real_t initRes = 0.01;
    bool converged = false;
    index_t iter = 0;

    gsStopwatch clock;
    clock.restart();
    // 0. pre-solve flow to get a stable flow
    gsNewton<real_t> newtonFlow(nsAssembler);
    newtonFlow.options().setInt("Verbosity",newton_verbosity::none);
    newtonFlow.options().setInt("Solver",linear_solver::LU);
    newtonFlow.solve();
    solutionFlow = newtonFlow.solution();
    nsAssembler.constructSolution(newtonFlow.solution(),newtonFlow.allFixedDofs(),velocity,pressure);

    while (!converged && iter < 50)
    {
        if (iter > 0)
        {
            // 1. compute ALE displacement
            aleAssembler.setFixedDofs(0,boundary::south,interfaceNow[0]-interfaceOld[0]);
            aleAssembler.setFixedDofs(1,boundary::north,interfaceNow[1]-interfaceOld[1]);
            aleAssembler.setFixedDofs(2,boundary::west,interfaceNow[2]-interfaceOld[2]);
            aleAssembler.assemble(ALE);
            gsSparseSolver<>::LU solverALE(aleAssembler.matrix());
            gsMatrix<> aleUpdateVector = solverALE.solve(aleAssembler.rhs());
            gsMultiPatch<> aleUpdate;
            aleAssembler.constructSolution(aleUpdateVector,aleUpdate);
            // 2. deform flow mesh
            // update ALE
            for (index_t p = 0; p < 3; ++p)
                ALE.patch(p).coefs() += aleUpdate.patch(p).coefs();
            // update flow mesh
            for (index_t p = 0; p < 3; ++p)
                nsAssembler.patches().patch(p+3).coefs() += aleUpdate.patch(p).coefs();
        }
        // 3. solve flow
        nsAssembler.assemble(velocity,pressure);
        gsSparseSolver<>::LU solverFlow(nsAssembler.matrix());
        solutionFlow += solverFlow.solve(nsAssembler.rhs());
        nsAssembler.constructSolution(solutionFlow,newtonFlow.allFixedDofs(),velocity,pressure);

        // 4. solve beam
        elAssembler.assemble(displacement);
        gsSparseSolver<>::LU solverBeam(elAssembler.matrix());
        solutionBeam += solverBeam.solve(elAssembler.rhs());
        elAssembler.constructSolution(solutionBeam,displacement);

        // 5. Aitken relaxation
        if (iter== 0)
        {
            interfaceNow.clear();
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::north)->coefs());
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::south)->coefs());
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::east)->coefs());
            omega = 1.;
        }
        else
        {
            interfaceNew.clear();
            interfaceNew.push_back(displacement.patch(0).boundary(boundary::north)->coefs());
            interfaceNew.push_back(displacement.patch(0).boundary(boundary::south)->coefs());
            interfaceNew.push_back(displacement.patch(0).boundary(boundary::east)->coefs());
            aitkenRelaxation(interfaceOld,interfaceNow,interfaceNew,omega);
        }

        // 7. plot
        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"fsi_steady_flow",collectionFlow,iter+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"fsi_steady_beam",collectionBeam,iter+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsPart,"fsi_steady_flow_part",collectionFlowPart,iter+1,numPlotPoints);
        }

        // 8. convergence
        real_t residual = computeResidual(interfaceOld,interfaceNow);
        if (iter == 0)
            initRes = residual;
        gsInfo << "Iter " << iter + 1 << ", intRes: " << residual << std::endl;
        if (residual < absTol || residual/initRes < relTol)
            converged = true;
        iter++;
    }
    gsInfo << "Solved in " << clock.stop() << "s.\n";

    // 5*. validation
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::east)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(4,boxSide(boundary::north)));
    bdrySides.push_back(std::pair<index_t,index_t>(5,boxSide(boundary::west)));
    gsMatrix<> force = nsAssembler.computeForce(velocity,pressure,bdrySides);
    gsInfo << "Drag: " << force.at(0) << std::endl;
    gsInfo << "Lift: " << force.at(1) << std::endl;
    gsMatrix<> point(2,1);
    point << 1.,0.5;
    gsInfo << "Displacement of the beam point A:\n" << displacement.patch(0).eval(point) << std::endl;

    gsInfo << "Plotting the output to the Paraview files \"fsi_steady flow.pvd\" and \"fsi_steady beam.pvd\"...\n";
    collectionFlow.save();
    collectionBeam.save();
    collectionFlowPart.save();
    return 0;
}
