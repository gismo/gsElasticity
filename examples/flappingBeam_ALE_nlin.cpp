/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
///
/// weak coupling
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Testing the steady fluid-structure interaction solver in 2D.\n";

    std::string filenameALE = ELAST_DATA_DIR"/flappingBeam_flowFull.xml";
    std::string filenameBeam = ELAST_DATA_DIR"/flappingBeam_beam.xml";
    index_t numUniRef = 3; // number of h-refinements for the fluid
    index_t numPlotPoints = 1000;
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t densitySolid = 1.0e3;
    real_t timeStep = 0.01;
    real_t timeSpan = 1;
    real_t meshPR = 0.3;
    real_t load = 1.;
    real_t stiff = 0.;
    index_t matLaw = 1;
    bool stop = true;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the steady fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications for the fluid",numUniRef);
    cmd.addInt("p","plot","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("m","mesh","Pois ratio for mesh",meshPR);
    cmd.addReal("l","load","Gravity",load);
    cmd.addReal("j","jac","Stiffening",stiff);
    cmd.addInt("n","law","Law",matLaw);
    cmd.addSwitch("b","break","Break the simulation",stop);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoALE; // this is a part of the flow geometry; we deform only this part to save memory and time
    gsReadFile<>(filenameALE, geoALE);
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

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

    // source function, rhs
    gsConstantFunction<> gBeam(0.,load*densitySolid,2);
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
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";
    // elasticity assembler: ALE
    gsElasticityAssembler<real_t> aleAssembler(geoALE,basisALE,bcInfoALE,gALE);
    aleAssembler.options().setReal("PoissonsRatio",meshPR);
    aleAssembler.options().setInt("MaterialLaw",matLaw);
    aleAssembler.options().setReal("LocalStiff",stiff);
    gsNewton<real_t> aleNewton(aleAssembler,gsMatrix<>::Zero(aleAssembler.numDofs(),1),aleAssembler.allFixedDofs());
    aleNewton.options().setInt("Verbosity",newton_verbosity::none);
    aleNewton.options().setInt("MaxIters",1);
    aleNewton.options().setInt("Solver",linear_solver::LDLT);
    gsInfo << "Initialized elasticity system for ALE with " << aleAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    gsMultiPatch<> displacement,ALE;
    // isogeometric fields (geometry + solution)
    gsField<> displacementField(geoBeam,displacement);
    gsField<> aleField(geoALE,ALE);
    gsPiecewiseFunction<> jacs;
    for (size_t p = 0; p < geoALE.nPatches(); ++p)
        jacs.addPiecePointer(new gsDetFunction<real_t>(geoALE,p));

    gsField<> jacField(geoALE,jacs,true);

    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;
    fieldsALE["Jac"] = &jacField;
    // paraview collection of time steps
    gsParaviewCollection collectionBeam("flappingBeam_mesh_beam");
    gsParaviewCollection collectionALE("flappingBeam_mesh_ALE");
    gsParaviewCollection collectionMesh("flappingBeam_mesh_mesh");

    gsProgressBar bar;
    gsStopwatch clock;

    std::ofstream file;
    file.open("flappingBeam_mesh_nonlin.txt");
    file << 0. << " " << 0. << std::endl;

    //=============================================//
             // Initial conditions //
    //=============================================//

    // plotting initial condition
    elAssembler.constructSolution(gsMatrix<>::Zero(elAssembler.numDofs(),1),displacement);
    aleAssembler.constructSolution(gsMatrix<>::Zero(aleAssembler.numDofs(),1),ALE);
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_mesh_beam",collectionBeam,0,numPlotPoints);
        gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_mesh_ALE",collectionALE,0,numPlotPoints);
        plotDeformation(geoALE,ALE,"flappingBeam_mesh_mesh",collectionMesh,0);
    }

    std::vector<gsMatrix<> > interfaceNow, interfaceNew;
    interfaceNow.push_back(ALE.patch(3).boundary(boundary::south)->coefs());
    interfaceNow.push_back(ALE.patch(4).boundary(boundary::north)->coefs());
    interfaceNow.push_back(ALE.patch(5).boundary(boundary::west)->coefs());
    interfaceNew = interfaceNow;

    elTimeSolver.initialize();

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    clock.restart();

    gsInfo << "Running the coupled simulation...\n";
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));
        //gsInfo << "Time " << i << std::endl;
        index_t bad = aleAssembler.checkSolution(ALE);

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

        // 4. ALE
        aleAssembler.setFixedDofs(3,boundary::south,interfaceNew[0]-interfaceNow[0]);
        aleAssembler.setFixedDofs(4,boundary::north,interfaceNew[1]-interfaceNow[1]);
        aleAssembler.setFixedDofs(5,boundary::west,interfaceNew[2]-interfaceNow[2]);
        aleNewton.reset();
        aleNewton.solve();
        // construct ALE
        aleAssembler.constructSolution(aleNewton.solution(),aleNewton.allFixedDofs(),ALE);

        interfaceNow = interfaceNew;

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flappingBeam_mesh_beam",collectionBeam,i+1,numPlotPoints);
            gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flappingBeam_mesh_ALE",collectionALE,i+1,numPlotPoints);
            plotDeformation(geoALE,ALE,"flappingBeam_mesh_mesh",collectionMesh,i+1);
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
        collectionMesh.save();
    }

    file.close();

    return 0;
}
