/// This is a simple numerical example of modeling active muscle behavior. It is based on the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
/// Muscle behavior is modeled with the incompressible nonlinear elasticity equations (pressure-displacement formulation)
/// with tendon and muscle materials for the passive part and a fiber-based active response part.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/src/gsMuscleAssembler.h>
#include <gsElasticity/src/gsMassAssembler.h>
#include <gsElasticity/src/gsElTimeIntegrator.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a simulation of active muscle behavior.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = gsElasticity_DATA"/muscleBeamMP.xml";
    real_t youngsModulusMuscle = 3.0e5; // shear modulus 1e5;
    real_t youngsModulusTendon = 3.0e6; // shear modulus 1e6;
    real_t poissonsRatioMuscle = 0.5;
    real_t poissonsRatioTendon = 0.5;
    real_t density = 9e2;
    real_t gravityAcc = 0.;
    real_t maxMuscleStress = 3.0e5; // maximum stress produced at the optimal fiber stretch
    real_t optFiberStretch = 1.3;   // optimal fiber stretch
    real_t deltaW = 0.3;            // shape parameter of the active reponse function
    real_t powerNu = 4.0;           // another shape parameter of the active reponse function
    // 1 = muscle, 0 = tendon
    gsFunctionExpr<> tendonMuscleSinglePatch("16*(1-x)^2*x^2",3);
    real_t prestress = 30;          // prestress applied at the ends of the muscle beam

    // direction of muscle fibers in the parametric domain
    gsVector<> fiberDirection(3);
    fiberDirection << 1.,0.,0.;
    // space discretization
    index_t numUniRefDirX = 2;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    bool subgridOrTaylorHood = false;
    // time integration
    real_t timeSpan = 2;
    real_t timeStep = 0.1;
    // output
    index_t numPlotPoints = 1000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is a simulation of active muscle behavior.");
    cmd.addInt("x","xrefine","Number of uniform refinement along the beam axis",numUniRefDirX);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation applications",numDegElev);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgridOrTaylorHood);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step, sec",timeStep);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("m","maxStress","Maximum stress produced at the optimal fiber stretch",maxMuscleStress);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << "Using " << (subgridOrTaylorHood ? "Taylor-Hood " : "subgrid ") << "mixed elements.\n";

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    geometry.computeTopology();

    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
    {
        basisDisplacement.degreeElevate();
        basisPressure.degreeElevate();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisDisplacement.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        for (index_t i = 0; i < numUniRefDirX; ++i)
        {
            basisDisplacement.basis(p).uniformRefine(1,1,0);
            basisPressure.basis(p).uniformRefine(1,1,0);
        }
    // additional displacement refinement for stable mixed FEM
    if (!subgridOrTaylorHood) // subgrid
        basisDisplacement.uniformRefine();
    else  // Taylor-Hood
        basisDisplacement.degreeElevate();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsConstantFunction<> prestressLeft(-1*prestress,0.,0.,3);
    gsConstantFunction<> prestressRight(prestress,0.,0.,3);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (size_t p = 0; p < geometry.nPatches(); ++p)
    {
        bcInfo.addCondition(p,boundary::west,condition_type::neumann,&prestressLeft);
        bcInfo.addCondition(p,boundary::east,condition_type::neumann,&prestressRight);
    }
    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,gravityAcc*density,3);

    gsPiecewiseFunction<> tendonMuscleDistribution;
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        tendonMuscleDistribution.addPiece(tendonMuscleSinglePatch);



    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // creating stiffness assembler
    gsMuscleAssembler<real_t> assembler(geometry,basisDisplacement,basisPressure,bcInfo,gravity,tendonMuscleDistribution,fiberDirection);
    assembler.options().setReal("MuscleYoungsModulus",youngsModulusMuscle);
    assembler.options().setReal("MusclePoissonsRatio",poissonsRatioMuscle);
    assembler.options().setReal("TendonYoungsModulus",youngsModulusTendon);
    assembler.options().setReal("TendonPoissonsRatio",poissonsRatioTendon);
    assembler.options().setReal("MaxMuscleStress",maxMuscleStress);
    assembler.options().setReal("OptFiberStretch",optFiberStretch);
    assembler.options().setReal("DeltaW",deltaW);
    assembler.options().setReal("PowerNu",powerNu);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(geometry,basisDisplacement,bcInfo,gravity);
    massAssembler.options().setReal("Density",density);

    // creating time integrator
    gsElTimeIntegrator<real_t> timeSolver(assembler,massAssembler);
    timeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    timeSolver.options().setInt("Verbosity",solver_verbosity::none);

    //=============================================//
            // Setting output & auxilary//
    //=============================================//

    // displacement field
    gsMultiPatch<> displacement, pressure;
    // stress field
    gsPiecewiseFunction<> stresses;
    // constructing an IGA field (geometry + solution)
    gsField<> dispField(geometry,displacement);
    gsField<> presField(geometry,pressure);
    gsField<> muscleTendonField(geometry,tendonMuscleDistribution,true);
    gsField<> stressField(assembler.patches(),stresses,true);
    std::map<std::string,const gsField<> *> fields;
    fields["Displacement"] = &dispField;
    fields["Pressure"] = &presField;
    fields["Muscle/tendon"] = &muscleTendonField;
    fields["von Mises"] = &stressField;

    gsProgressBar bar;
    gsStopwatch totalClock;

    //=============================================//
                   // Initial conditions //
    //=============================================//

    // set initial conditions
    timeSolver.setDisplacementVector(gsMatrix<>::Zero(massAssembler.numDofs(),1));
    timeSolver.setVelocityVector(gsMatrix<>::Zero(massAssembler.numDofs(),1));

    timeSolver.constructSolution(displacement,pressure);
    assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

    // plotting initial displacement
    gsParaviewCollection collection("muscleBeam");
    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields,"muscleBeam",collection,0,numPlotPoints);

    //=============================================//
                  // Solving //
    //=============================================//

    gsInfo << "Running the simulation...\n";
    totalClock.restart();
    for (index_t i = 0; i < index_t(timeSpan/timeStep); ++i)
    {
        bar.display(i+1,index_t(timeSpan/timeStep));

        assembler.options().setReal("Alpha",(1-cos(EIGEN_PI*(i+1)*timeStep))/2);
        timeSolver.makeTimeStep(timeStep);
        // construct solution; timeSolver already knows the new Dirichlet BC
        timeSolver.constructSolution(displacement,pressure);
        assembler.constructCauchyStresses(displacement,pressure,stresses,stress_components::von_mises);

        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields,"muscleBeam",collection,i+1,numPlotPoints);
    }

    //=============================================//
                // Final touches //
    //=============================================//

    gsInfo << "Simulation time: " + secToHMS(totalClock.stop()) + "\n";

    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"muscleBeam.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
