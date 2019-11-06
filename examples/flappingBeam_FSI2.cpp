/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
///
/// initial version, september 2019
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
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

    for (index_t i = 0; i < interfaceOld.size(); ++i)
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
    gsInfo << "Testing the steady fluid-structure interaction solver in 2D.\n";

    std::string filenameFlow = ELAST_DATA_DIR"/fsi_flow_around_cylinder.xml";
    std::string filenameFlowPart = ELAST_DATA_DIR"/fsi_flow_around_cylinder_segment.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/fsi_beam_around_cylinder.xml";
    index_t numUniRefFlow = 2; // number of h-refinements for the fluid
    index_t numKRefFlow = 0; // number of k-refinements for the fluid
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numUniRefBeam = 2; // number of h-refinements for the beam and the ALE mapping
    index_t numKRefBeam = 0; // number of k-refinements for the beam and the ALE mapping
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    bool subgrid = true;
    bool supg = false;
    real_t densityFluid = 1.0e3;
    real_t densitySolid = 1.0e4;
    real_t absTol = 1e-6;
    real_t relTol = 1e-6;
    real_t timeStep = 0.001;
    real_t timeSpan = 1.;
    real_t warmUpTimeSpan = 3.;
    real_t warmUpTimeStep = 0.1;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the steady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRefFlow);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("b","beamrefine","Number of uniform refinement applications for the beam and ALE",numUniRefBeam);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    gsMultiPatch<> geoPart; // this is a part of the flow geometry; we deform only this part to save memory and time
    gsReadFile<>(filenameFlowPart, geoPart);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

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
    basisDisplacement.uniformRefine();
    geoPart.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisALE(geoPart);

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

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

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,g);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,g);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler);
    nsTimeSolver.options().setInt("Scheme",time_integration_NS::theta_scheme_linear);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: flow mesh
    gsElasticityAssembler<real_t> aleAssembler(geoPart,basisALE,bcInfoALE,g);
    aleAssembler.options().setReal("PoissonsRatio",0.4);
    aleAssembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsInfo << "Initialized elasticity system for ALE with " << aleAssembler.numDofs() << " dofs.\n";

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
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flappingBeam_FSI2_flow");
    gsParaviewCollection collectionBeam("flappingBeam_FSI2_beam");
    gsParaviewCollection collectionALE("flappingBeam_FSI2_ALE");

    gsProgressBar bar;
    gsStopwatch clock;

    //=============================================//
                   // Warming up //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.getFixedDofs(0,boundary::west,inflowDDoFs);

    // set all Dirichlet DoFs to zero
    nsAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());
    nsTimeSolver.initialize();

    clock.restart();
    gsInfo << "Running the simulation with a coarse time step to compute an initial solution...\n";
    for (index_t i = 0; i < index_t(warmUpTimeSpan/warmUpTimeStep); ++i)
    {
        bar.display(i+1,index_t(warmUpTimeSpan/warmUpTimeStep));
        if ((i+1)*warmUpTimeStep < 2./3.*warmUpTimeSpan)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(M_PI*warmUpTimeStep*(i+1)/2./warmUpTimeSpan*3.))/2);
        else
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs);
        nsTimeSolver.makeTimeStep(warmUpTimeStep);
    }
    gsInfo << "Complete in " << clock.stop() << "s.\n";

    //=============================================//
             // Initial conditions //
    //=============================================//



    gsMatrix<> solutionFlow = nsTimeSolver.solutionVector();
    gsMatrix<> solutionBeam = gsMatrix<>::Zero(elAssembler.numDofs(),1);
    gsMatrix<> solutionALE = gsMatrix<>::Zero(aleAssembler.numDofs(),1);

    // plotting initial condition
    nsAssembler.constructSolution(solutionFlow,velocity,pressure);
    elAssembler.constructSolution(solutionBeam,displacement);
    aleAssembler.constructSolution(solutionALE,ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"fsi_unsteady_flow",collectionFlow,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"fsi_unsteady_beam",collectionBeam,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsPart,"fsi_unsteady_flow_part",collectionFlowPart,0,numPlotPoints);
    }

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


    //=============================================//
             // Coupled simulation //
    //=============================================//

    // container for interface displacement
    std::vector<gsMatrix<> > interface;
    interface.push_back(ALE.patch(0).boundary(boundary::south)->coefs());
    interface.push_back(ALE.patch(1).boundary(boundary::north)->coefs());
    interface.push_back(ALE.patch(2).boundary(boundary::west)->coefs());
    gsMultiPatch<> updateALE;
    aleAssembler.constructSolution(solutionALE,updateALE);
    std::vector<std::pair<index_t,index_t> > alePatches;
    alePatches.push_back(std::pair<index_t,index_t>(3,0));
    alePatches.push_back(std::pair<index_t,index_t>(4,1));
    alePatches.push_back(std::pair<index_t,index_t>(5,2));
    nsAssembler.setALEvelocity(updateALE,alePatches);


    gsMatrix<> solutionBeamVelocity = gsMatrix<>::Zero(elAssembler.numDofs(),1);
    gsMultiPatch<> trueVelocity;
    nsAssembler.constructSolution(solutionFlow,trueVelocity);
    gsMatrix<> tempFlowSOl;
    // time loop
    clock.restart();
    for (index_t i = 0; i < 1; ++i)
    {
        gsInfo << "=========================================================TIME STEP " << i+1 << "/" << index_t(timeSpan/timeStep) << std::endl;
        std::vector<gsMatrix<> > interfaceOld, interfaceNow, interfaceNew;
        interfaceOld.push_back(interface[0]);
        interfaceOld.push_back(interface[1]);
        interfaceOld.push_back(interface[2]);
        // Aitken relaxation factor
        real_t omega;
        // convergence
        real_t initRes;
        bool converged = false;
        index_t iter = 0;

        // fsi coupling loop
        while (!converged && iter < 3 )
        {
            if (iter > 0)
            {
                /*for (index_t p = 0; p < 3; ++p)
                {
                    nsAssembler.patches().patch(p+3).coefs() -= updateALE.patch(p).coefs()*timeStep;
                    ALE.patch(p).coefs() -= updateALE.patch(p).coefs()*timeStep;
                }

                // 1. compute ALE displacement
                aleAssembler.setDirichletDofs(0,boundary::south,interfaceNow[0]-interface[0]);
                aleAssembler.setDirichletDofs(1,boundary::north,interfaceNow[1]-interface[1]);
                aleAssembler.setDirichletDofs(2,boundary::west,interfaceNow[2]-interface[2]);
                aleAssembler.assemble(ALE);
                gsSparseSolver<>::LU solverALE(aleAssembler.matrix());

                gsMatrix<> vectorUpdateALE = solverALE.solve(aleAssembler.rhs());




                aleAssembler.constructSolution(vectorUpdateALE,updateALE);
                for (index_t p = 0; p < 3; ++p)
                    updateALE.patch(p).coefs() /= timeStep;
            }

            // 2. solve flow
            // deform flow mesh
            for (index_t p = 0; p < 3; ++p)
            {
                nsAssembler.patches().patch(p+3).coefs() += updateALE.patch(p).coefs()*timeStep;
                ALE.patch(p).coefs() += updateALE.patch(p).coefs()*timeStep;
            }

           //nsAssembler.setDirichletDofs(3,boundary::south,
             //                            updateALE.patch(0).boundary(boundary::south)->coefs());
            //nsAssembler.setDirichletDofs(4,boundary::north,
              //                           updateALE.patch(1).boundary(boundary::north)->coefs());
            //nsAssembler.setDirichletDofs(5,boundary::west,
              //                           updateALE.patch(2).boundary(boundary::west)->coefs());

            //nsTimeSolver.setInitialSolution(solutionFlow);

            //nsTimeSolver.makeTimeStep(timeStep);
       tempFlowSOl =  nsTimeSolver.oseenFSI(trueVelocity,timeStep,solutionFlow);

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

            // 4. solve beam
            gsInfo << "beam\n";
            elTimeSolver.setInitialDisplacement(solutionBeam);
            elTimeSolver.setInitialVelocity(solutionBeamVelocity);
            elTimeSolver.initialize();


            // set new forces AFTER the initialization where the old residual is assembled
            nsAssembler.constructSolution(tempFlowSOl,velocity,pressure);


            elTimeSolver.makeTimeStep(timeStep);
            elAssembler.constructSolution(elTimeSolver.displacementVector(),displacement);
            gsMatrix<> point(2,1);
            point << 1.,0.5;
            gsInfo << "Displacement of the beam point A:\n" << displacement.patch(0).eval(point) << std::endl;
            gsMatrix<> pointA(2,1);
            pointA << 0.,0.5;
            gsInfo << "ALE A:\n" << ALE.patch(2).eval(pointA) << std::endl;

            if (iter > 0)
            {
                interfaceOld.clear();
                interfaceOld.push_back(interfaceNow[0]);
                interfaceOld.push_back(interfaceNow[1]);
                interfaceOld.push_back(interfaceNow[2]);

            }
            interfaceNow.clear();
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::north)->coefs());
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::south)->coefs());
            interfaceNow.push_back(displacement.patch(0).boundary(boundary::east)->coefs());

             // 5. Aitken relaxation
            if (iter == 0)
            {
                interfaceNew.clear();
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

            // 8. convergence
            real_t residual = computeResidual(interfaceOld,interfaceNow);
            if (iter == 0)
                initRes = residual;
            gsInfo << "Iter " << iter + 1 << ", absRes: " << residual
                   << ", relRes: " << residual/initRes << std::endl;
            if (residual < absTol || residual/initRes < relTol)
                converged = true;

            if (numPlotPoints > 0)
            {
                gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"fsi_insteady_flow",collectionFlow,iter+1,numPlotPoints);
                gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"fsi_unsteady_beam",collectionBeam,iter+1,numPlotPoints);
                gsWriteParaviewMultiPhysicsTimeStep(fieldsPart,"fsi_unsteady_flow_part",collectionFlowPart,iter+1,numPlotPoints);
            }

            ++iter;

        }



        for (index_t i = 0; i < 3; ++i)
            interface[i] = interfaceNow[i];

        solutionBeam = elTimeSolver.displacementVector();
        solutionBeamVelocity = elTimeSolver.velocityVector();
        nsAssembler.constructSolution(nsTimeSolver.solutionVector(),trueVelocity);
        solutionFlow = tempFlowSOl;

    }


    gsInfo << "Solved in " << clock.stop() << "s.\n";

    gsInfo << "The output is plotted to the Paraview files \"fsi_unsteady flow.pvd\" and \"fsi_unsteady beam.pvd\"...\n";
    collectionFlow.save();
    collectionBeam.save();
    collectionFlowPart.save();

    return 0;*/
}
