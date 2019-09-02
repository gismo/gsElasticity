/// This is an example of using the Navier-Stokes solver on a 2D multi-patch geometry
#include <gismo.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNewton.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){
    gsInfo << "Testing the Navier-Stokes solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/flow_around_cylinder.xml";
    index_t numUniRef = 3; // number of h-refinements
    index_t numKRef = 0; // number of k-refinements
    index_t numBLRef = 1; // number of additional boundary layer refinements for the fluid
    index_t numPlotPoints = 10000;
    real_t viscosity = 0.001;
    real_t maxInflow = 0.3;
    bool subgrid = false;
    index_t iters = 20;
    bool supg = false;

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the Stokes solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addInt("k","krefine","Number of k refinement applications",numKRef);
    cmd.addInt("l","blayer","Number of additional boundary layer refinements for the fluid",numBLRef);
    cmd.addInt("s","sample","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("v","viscosity","Viscosity of the fluid",viscosity);
    cmd.addReal("f","inflow","Maximum inflow velocity",maxInflow);
    cmd.addSwitch("e","element","True - subgrid, false - TH",subgrid);
    cmd.addInt("i","iters","Max number of Newton's iterations",iters);
    cmd.addSwitch("g","supg","Do NOT use SUPG stabilization",supg);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);
    // inflow velocity profile
    gsFunctionExpr<> inflow(util::to_string(maxInflow) + "*4*y*(0.41-y)/0.41^2",2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,&inflow,0);
    bcInfo.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfo.addCondition(0,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfo.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
    }

    //=============================================//
                  // Assembly //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    // creating bases
    gsMultiBasis<> basisVelocity(geometry);
    gsMultiBasis<> basisPressure(geometry);
    for (index_t i = 0; i < numKRef; ++i)
    {
        basisVelocity.degreeElevate();
        basisPressure.degreeElevate();
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    for (index_t i = 0; i < numUniRef; ++i)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }
    // additional refinement of the boundary layer around the cylinder
    gsMatrix<> box(2,2);
    box << 0.,0.,0.,0.2;
    for (index_t i = 0; i < numBLRef; ++i)
        for (index_t p = 0; p < 4; ++p)
        {
            basisVelocity.refine(p,box);
            basisPressure.refine(p,box);
        }
    // additional velocity refinement for stable mixed FEM
    if (subgrid)
        basisVelocity.uniformRefine();
    else
        basisVelocity.degreeElevate();

    // creating assembler
    gsNsAssembler<real_t> assembler(geometry,basisVelocity,basisPressure,bcInfo,g);
    assembler.options().setReal("Viscosity",viscosity);
    assembler.options().setInt("DirichletValues",dirichlet::interpolation);
    assembler.options().setSwitch("SUPG",supg);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving Newton//
    //=============================================//

    // setting Newton's method
    gsNewton<real_t> newton(assembler);
    newton.options().setInt("Verbosity",newton_verbosity::all);
    newton.options().setInt("MaxIters",iters);
    newton.options().setInt("Solver",linear_solver::LU);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    newton.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                  // Output //
    //=============================================//

    // constructing solution as an IGA function
    gsMultiPatch<> velocity, pressure;
    assembler.constructSolution(newton.solution(),velocity,pressure);
    // constructuin isogeometric field (geometry + solution)
    gsField<> velocityField(assembler.patches(),velocity);
    gsField<> pressureField(assembler.patches(),pressure);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string,const gsField<> *> fields;
    fields["Velocity"] = &velocityField;
    fields["Pressure"] = &pressureField;
    gsInfo << "Plotting the output to the Paraview file \"NS_around_cylinder.pvd\"...\n";
    gsWriteParaviewMultiPhysics(fields,"NS_aroundCylinder",numPlotPoints);

    //=============================================//
                  // Validation //
    //=============================================//

    // computing force acting on the surface of the cylinder
    std::vector<std::pair<index_t, boxSide> > bdrySides;
    bdrySides.push_back(std::pair<index_t,index_t>(0,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(1,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(2,boxSide(boundary::south)));
    bdrySides.push_back(std::pair<index_t,index_t>(3,boxSide(boundary::south)));
    gsMatrix<> force = assembler.computeForce(velocity,pressure,bdrySides);
    real_t L = 0.1; // characteristic lenght
    real_t U_mean = maxInflow * 2./3.; // mean velocity
    gsInfo << "Drag coefficient: " << 2.*force.at(0)/L/pow(U_mean,2) << std::endl;
    gsInfo << "Lift coefficient: " << 2.*force.at(1)/L/pow(U_mean,2) << std::endl;

    // evaluating pressure difference at the far front and the far rear points of the cylinder
    gsMatrix<> point(2,1);
    point << 0.5, 0;                                                                   // this info is hard-cored in the geometry
    gsInfo << "Pressure difference: " << pressure.patch(0).eval(point)(0,0) -          // far front point
                                         pressure.patch(2).eval(point)(0,0) << "Pa\n"; // far rear point

    return 0;
}
