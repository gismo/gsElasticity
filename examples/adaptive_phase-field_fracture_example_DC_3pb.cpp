/** @file adaptive_phase-field_fracture_example_DC_3pb.cpp

    @brief Three-point bending variant of adaptive_phase-field_fracture_example_DC.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

    -----------------------------------------------------------------------
    Three-point bending (3PB) test on a rectangular beam.

    The supports and the load are NOT located on patch boundaries (sides),
    so they cannot be expressed with gsBoundaryConditions, which only
    understands whole patch sides/corners/interfaces. Instead:

      - The beam is kept as a SINGLE patch (so the existing adaptive THB
        refinement machinery, which assumes mb.nBases()==1, keeps working
        unmodified).
      - The left support   (x = supportLeftRatio  * L, bottom side, y=0)
      - the right support   (x = supportRightRatio * L, bottom side, y=0)
      - the loading point   (x = loadRatio         * L, top side,    y=H)
        are located by finding, among the basis functions whose support
        lies on the corresponding side, those whose anchor (Greville point)
        maps to the physical location closest to the target x -- in 2D this
        is a single dof (a side is a curve); in 3D it is generally a whole
        "row" of dofs spanning the full depth (z) of that side, matching a
        cylindrical roller support/loading nose spanning the specimen width.
      - Those degrees of freedom are then constrained directly by algebraic
        row/column elimination on the assembled linear system (a manual,
        exact Dirichlet-style elimination), applied every time the
        elasticity operator is (re-)assembled.

    See the accompanying design notes for why gsBoundaryConditions and
    gsDofMapper::eliminateDof() (the natural, but here inapplicable,
    built-in mechanisms) cannot be used without further changes to
    gsSolidAssembler, and why the row/column-surgery approach below needs
    no modification whatsoever to G+Smo or gsElasticity.

    To run the script in 2D:
    ./bin/adaptive_phase-field_fracture_example_DC_3pb -I <inputDir> --plot
*/

//! [Include namespace]
#include <gismo.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsMaterialEval.h>
#include <gsElasticity/gsSolidAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPSOR.h>
#include <gsUtils/gsStopwatch.h>
#include <gsHSplines/gsHElementHelper.h>
#include <gsHSplines/gsHElementMarker.h>

using namespace gismo;
//! [Include namespace]

#define PRINT(w) std::setw(w)<<std::left

std::vector<real_t> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const real_t         & lowerBound=0.0,
                                    const real_t         & upperBound=1.0);

template <short_t dim,class T>
T refineMesh(      gsMultiBasis<T> & basis,
             const std::vector<T>  & vals,
             const gsOptionList    & options = gsOptionList());

template <class T>
struct times
{
    T elAssemblyTime;
    T elSolverTime;
    T pfAssemblyTime;
    T pfSolverTime;
    T projectionTime;

    void reset()
    {
        elAssemblyTime = 0;
        elSolverTime = 0;
        pfAssemblyTime = 0;
        pfSolverTime = 0;
        projectionTime = 0;
    }
};

//////////////////////////////////////////////////////////////////////////
// 3-POINT-BENDING POINT-CONSTRAINT HELPERS ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

/// A single scalar constraint: fix component \a component of the displacement
/// at local basis function \a localDof of \a patch to \a value.
template <class T>
struct gsPointConstraint
{
    index_t     patch;
    index_t     localDof;
    index_t     component;
    T           value;
    gsVector<T> location; // physical location, kept for logging/plotting only
};

/// Physical bounding box of patch 0, assuming a straight-edged (bilinear or
/// higher order but axis-aligned) rectangular beam, so that the extremal
/// control points coincide with the extremal physical coordinates.
/// For curved/skewed geometries, replace this by evaluating the geometry
/// at the parametric corners instead of using the control net.
template <class T>
void getBoundingBox(const gsMultiPatch<T> & mp, gsVector<T> & lower, gsVector<T> & upper)
{
    const gsMatrix<T> & coefs = mp.patch(0).coefs();
    lower = coefs.colwise().minCoeff().transpose();
    upper = coefs.colwise().maxCoeff().transpose();
}

/** \brief Locate ALL local dofs (basis functions) on a given patch side whose
 *  anchor (Greville abscissa) maps to a physical location with x-coordinate
 *  closest to \a targetX -- i.e. the "control-point row"/knot line at x =
 *  targetX. In 2D this row degenerates to a single dof (a side is a curve);
 *  in 3D it is generally a whole line of dofs spanning the full depth (z)
 *  of that side, which is exactly what a cylindrical roller support/loading
 *  nose spanning the width of a 3PB specimen constrains.
 *
 * Implementation notes (this answers points 1-3 of the investigation):
 *  - gsBasis::boundary(side) returns the LOCAL indices of all basis
 *    functions with support on that side (already used in the original
 *    example to compute reaction forces).
 *  - gsBasis::anchors() returns, for every local dof (same indexing as
 *    boundary()/active()), its anchor point in parameter space.
 *  - gsGeometry::eval_into() maps that parameter point to a physical
 *    location, so that we can compare against a physical target
 *    coordinate regardless of the parametrization used.
 *  - This works unmodified for adaptively refined THB-spline bases
 *    (gsHTensorBasis implements both boundary() and anchors_into()).
 *
 * \a tol is an absolute physical-distance tolerance (same units as the
 * geometry) used to decide which dofs belong to the same "row" as the
 * closest one found; it should be well below half the local element size.
 */
template <class T>
std::vector<index_t> locateBoundaryLineAtX(const gsMultiBasis<T> & mb,
                                            const gsMultiPatch<T> & mp,
                                            index_t patch,
                                            boxSide side,
                                            T targetX,
                                            gsVector<T> & repLocation,
                                            T tol = (T)1e-6)
{
    gsMatrix<index_t> bnd     = mb.basis(patch).boundary(side);
    gsMatrix<T>        anchors = mb.basis(patch).anchors();

    std::vector<T> dist(bnd.rows());
    T bestDist = std::numeric_limits<T>::max();
    gsMatrix<T> pt, phys;
    for (index_t k = 0; k!=bnd.rows(); ++k)
    {
        pt = anchors.col(bnd(k,0));
        mp.patch(patch).eval_into(pt,phys);
        dist[k] = math::abs(phys(0,0)-targetX);
        if (dist[k]<bestDist)
        {
            bestDist     = dist[k];
            repLocation  = phys.col(0);
        }
    }
    GISMO_ASSERT(bestDist<std::numeric_limits<T>::max(),
                 "No boundary dof found on side "<<side<<" of patch "<<patch);

    std::vector<index_t> result;
    for (index_t k = 0; k!=bnd.rows(); ++k)
        if (dist[k] <= bestDist+tol)
            result.push_back(bnd(k,0));
    return result;
}

/// From a line of local dofs (as returned by locateBoundaryLineAtX), pick the
/// single one physically closest to \a targetCoord in component \a comp.
/// Used to reduce a line support/load to a single "pin" point (e.g. to break
/// the remaining rigid-body modes in 3D, or to model a true point indenter
/// instead of a line-support/roller).
template <class T>
index_t pickDofNearCoord(const std::vector<index_t> & line,
                         const gsMultiPatch<T> & mp, index_t patch,
                         const gsMatrix<T> & anchors,
                         T targetCoord, index_t comp)
{
    index_t best = -1;
    T bestDist = std::numeric_limits<T>::max();
    gsMatrix<T> pt, phys;
    for (const index_t dof : line)
    {
        pt = anchors.col(dof);
        mp.patch(patch).eval_into(pt,phys);
        const T d = math::abs(phys(comp,0)-targetCoord);
        if (d<bestDist) { bestDist = d; best = dof; }
    }
    GISMO_ASSERT(best>=0, "Empty dof line passed to pickDofNearCoord");
    return best;
}

/** \brief Build the two support constraints of a 3-point-bending test.
 *
 * Left support (x = xLeft) is modelled as a pin: u_y=0 along the WHOLE
 * support line (roller behaviour, all z), PLUS u_x=0 (and, in 3D, u_z=0)
 * at a single point on that line (its center in z) to remove the
 * remaining horizontal/out-of-plane rigid-body modes -- the standard,
 * minimal-constraint choice for 3-point-bending FE models.
 * Right support (x = xRight) is a pure roller: u_y=0 along the whole line.
 */
template <class T>
std::vector<gsPointConstraint<T>> findSupportDofs(const gsMultiBasis<T> & mb,
                                                    const gsMultiPatch<T> & mp,
                                                    T xLeft, T xRight,
                                                    T zMid, index_t dim)
{
    std::vector<gsPointConstraint<T>> result;
    const gsMatrix<T> anchors = mb.basis(0).anchors();
    gsVector<T> loc;

    const std::vector<index_t> lineL = locateBoundaryLineAtX(mb,mp,0,boundary::south,xLeft,loc);
    for (const index_t dof : lineL)
        result.push_back({0,dof,1,(T)0.0,loc}); // roller behaviour: u_y=0 everywhere on the line

    const index_t pinDof = (dim==3) ? pickDofNearCoord(lineL,mp,0,anchors,zMid,2) : lineL.front();
    result.push_back({0,pinDof,0,(T)0.0,loc}); // pin: u_x=0 at a single point only
    if (dim==3)
        result.push_back({0,pinDof,2,(T)0.0,loc}); // pin: u_z=0 at the same point (removes z rigid-body mode)

    const std::vector<index_t> lineR = locateBoundaryLineAtX(mb,mp,0,boundary::south,xRight,loc);
    for (const index_t dof : lineR)
        result.push_back({0,dof,1,(T)0.0,loc}); // roller: u_y=0 only

    return result;
}

/// Build the loading-line constraint: prescribe the vertical displacement
/// at x = xLoad on the TOP side to -uy (downward, displacement-controlled),
/// along the whole width (z) of the specimen -- i.e. a line-load, matching a
/// loading nose/roller spanning the specimen width. Set \a asPoint=true to
/// instead apply it only at the single dof closest to z=zMid (a true point
/// indenter).
template <class T>
std::vector<gsPointConstraint<T>> findLoadDofs(const gsMultiBasis<T> & mb,
                                                const gsMultiPatch<T> & mp,
                                                T xLoad, T uy,
                                                T zMid, index_t dim,
                                                bool asPoint = false)
{
    std::vector<gsPointConstraint<T>> result;
    gsVector<T> loc;
    const std::vector<index_t> line = locateBoundaryLineAtX(mb,mp,0,boundary::north,xLoad,loc);

    if (dim==3 && asPoint)
    {
        const gsMatrix<T> anchors = mb.basis(0).anchors();
        const index_t dof = pickDofNearCoord(line,mp,0,anchors,zMid,2);
        result.push_back({0,dof,1,-uy,loc});
    }
    else
    {
        for (const index_t dof : line)
            result.push_back({0,dof,1,-uy,loc});
    }
    return result;
}

/** \brief Build a gsDofMapper whose numbering is guaranteed to match the
 *  internal numbering used by a gsSolidAssembler constructed with \a bc,
 *  PROVIDED that \a bc only contains "Dirichlet" (side) conditions (no
 *  corners/collapsed/coupled conditions, no multi-patch interfaces).
 *
 * This mirrors exactly the "Dirichlet" branch of gsFeSpace<T>::setup()
 * (src/gsExpressions/gsFeSpace.h): gsDofMapper(mb,dim) followed by
 * markBoundary() for every Dirichlet side, then finalize(). Since our
 * point constraints are NOT part of \a bc, they remain free dofs in this
 * mapper, exactly as they remain free dofs in the assembler's own mapper;
 * we eliminate them ourselves afterwards, algebraically, in
 * applyPointDirichlet().
 */
template <class T>
gsDofMapper buildCompanionMapper(const gsMultiBasis<T> & mb,
                                  const gsBoundaryConditions<T> & bc,
                                  index_t dim)
{
    gsDofMapper mapper(mb,dim);
    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin("Dirichlet"); it!=bc.end("Dirichlet"); ++it)
    {
        bnd = mb.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, it->unkComponent());
    }
    mapper.finalize();
    return mapper;
}

/** \brief Impose a set of point-wise Dirichlet constraints on an already
 *  assembled linear system (K,F) by exact algebraic row/column elimination.
 *
 * For every constrained global dof d with prescribed value v:
 *   1. its column's contribution is moved to the right-hand side of every
 *      OTHER (still free) equation:  F(i) -= K(i,d)*v  for i != d;
 *   2. row d and column d of K are zeroed, K(d,d) is set to 1;
 *   3. F(d) is set to v.
 *
 * The resulting system, solved for the full vector u, satisfies u(d)=v
 * exactly (up to solver tolerance) -- this is standard Dirichlet
 * elimination performed as a post-assembly step, entirely in user code.
 * It requires no change to gsExprAssembler, gsDofMapper or gsSolidAssembler.
 *
 * This is the manual counterpart of what gsExprAssembler does internally
 * for ordinary boundary Dirichlet dofs (see gsExpressions/gsFeSpace.h,
 * the `push` routine referenced in the design notes): there, dofs marked
 * via gsDofMapper::eliminateDof() are removed from the system and their
 * contribution is moved to the rhs of the free rows. Here we cannot reach
 * that mechanism for interior dofs without changing gsSolidAssembler (see
 * design notes), so we reproduce its effect directly on the assembled
 * (full-size) matrix and vector.
 */
template <class T>
void applyPointDirichlet(gsSparseMatrix<T> & K,
                          gsMatrix<T>       & F,
                          const gsDofMapper & mapper,
                          const std::vector<gsPointConstraint<T>> & dofs)
{
    std::vector<index_t> gidx;
    std::vector<T>       gval;
    gidx.reserve(dofs.size());
    gval.reserve(dofs.size());
    for (const gsPointConstraint<T> & pc : dofs)
    {
        gidx.push_back(mapper.index(pc.localDof, pc.patch, pc.component));
        gval.push_back(pc.value);
    }

    // 1) Move known column contributions to the rhs (using the still-intact K)
    for (size_t c = 0; c!=gidx.size(); ++c)
    {
        const index_t d = gidx[c];
        const T       v = gval[c];
        for (typename gsSparseMatrix<T>::InnerIterator it(K,d); it; ++it)
            if (it.row()!=d)
                F(it.row(),0) -= it.value()*v;
    }

    // 2) Zero the rows/columns of all constrained dofs, set diagonal to 1
    std::vector<bool> isConstrained(K.rows(),false);
    for (index_t d : gidx) isConstrained[d] = true;

    for (index_t k = 0; k<K.outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(K,k); it; ++it)
            if (isConstrained[it.row()] || isConstrained[it.col()])
                it.valueRef() = (it.row()==it.col() && isConstrained[it.row()]) ? T(1) : T(0);

    // 3) Prescribe the values on the rhs
    for (size_t c = 0; c!=gidx.size(); ++c)
        F(gidx[c],0) = gval[c];
}

//////////////////////////////////////////////////////////////////////////

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & solverParameters,
           gsOptionList & mesherOptions,
           gsMultiPatch<T> & mp,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           std::string & outputdir);

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numHRef = 0;
    index_t numElev = 0;
    std::string output;
    std::string parInput;
    std::string geoInput;
    std::string damageInput;
    std::string inputDir;

    gsCmdLine cmd("Tutorial on solving a three-point-bending phase-field fracture problem.");
    cmd.addInt("e", "numElev","Number of degree elevation steps to perform before solving",numElev);
    cmd.addInt("r", "numHRef","Number of Uniform h-refinement loops", numHRef);
    cmd.addInt("p", "plotmod","Modulo for plotting", plotmod);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("g", "geometry", "Geometry file", geoInput);
    cmd.addString("d", "damage", "Damage file", damageInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    inputDir = inputDir + gsFileManager::getNativePathSeparator();
    std::string parInputPath = (parInput.empty() ? inputDir + "parameters.xml" : parInput);
    std::string geoInputPath = (geoInput.empty() ? inputDir + "geometry.xml" : geoInput);
    std::string damageInputPath = (damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    GISMO_ASSERT(gsFileManager::fileExists(parInputPath), "Input parameter file "<<parInputPath<<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(geoInputPath), "Input geometry file "<<geoInputPath<<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(damageInputPath), "Input damage file "<<damageInputPath<<" not found.");
    gsInfo << "Input parameter file "<< parInputPath <<"\n";
    gsInfo << "Input geometry file "<< geoInputPath <<"\n";
    gsInfo << "Input damage file "<< damageInputPath <<"\n";

    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    gsFileData<> fd_geo(geoInput.empty() ? inputDir + "geometry.xml" : geoInput);
    gsMultiPatch<> mp_ini;
    fd_geo.getFirst(mp_ini);
    if (numElev > 0)
        mp_ini.degreeIncrease(numElev);
    for (index_t i = 0; i<numHRef; ++i)
        mp_ini.uniformRefine(1);
    if (plot) gsWriteParaview(mp_ini,outputdir+"mp",10,true);

    gsFileData<> fd_damage(damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    gsMultiPatch<> damage;
    fd_damage.getFirst(damage);
    if (plot) gsWriteParaview(mp_ini,damage,outputdir+"initial_damage",100000);

    gsFileData<> fd_pars(parInput.empty() ? inputDir + "parameters.xml" : parInput);
    GISMO_ASSERT(fd_pars.hasLabel("material"), "Material parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("control"), "Displacement-Control parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("meshing"), "Adaptive meshing parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_u"), "Displacement boundary conditions not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_d"), "Phase-field boundary conditions not found in the input file.");

    //// Material parameters
    gsOptionList materialParameters;
    fd_pars.getLabel("material", materialParameters);

    //// Boundary control parameters
    gsOptionList controlParameters;
    fd_pars.getLabel("control", controlParameters);

    //// Solver parameters
    gsOptionList solverParameters;
    if (fd_pars.hasLabel("solver"))
        fd_pars.getLabel("solver", solverParameters);

    //// Boundary conditions
    // NOTE: for a pure 3-point-bending setup, BCs_u is typically EMPTY:
    // none of the supports/load sit on a whole patch side. It is still
    // read (and passed through) so that additional, ordinary side
    // conditions (e.g. out-of-plane symmetry in a 3D model) remain
    // possible without any change to this example.
    gsBoundaryConditions<> bc_u;
    fd_pars.getLabel("BCs_u", bc_u);

    //// Boundary conditions
    gsBoundaryConditions<> bc_d;
    fd_pars.getLabel("BCs_d", bc_d);

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Call the dimensional solver
    ///////////////////////////////////////////////////////////////////////////////////////
    switch (mp_ini.domainDim())
    {
        case 2:
            solve<2>(materialParameters,controlParameters,solverParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,solverParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        default:
            GISMO_ERROR("Invalid domain dimension");
    }

    return 0;
} // end main

template<short_t dim, class T>
std::vector<T> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const T         & lowerBound,
                                    const T         & upperBound)
{
    GISMO_ASSERT(basis.nBases() == 1, "Labeling is only implemented for single basis meshes");
    typename gsBasis<T>::domainIter domEnd = basis.basis(0).domain()->endAll();
    std::vector<T> labels(basis.basis(0).numElements());
    gsVector<index_t,dim> np;
    np.setConstant(2);
    gsLobattoRule<T> rule(np); // equivalent to using gsPointGrid with np = 2
    for (typename gsBasis<T>::domainIter domIt  = basis.basis(0).domain()->beginAll(); domIt<domEnd; ++domIt)
    {
        gsMatrix<T> nodes;
        gsMatrix<T> vals;
        gsVector<T> weights;

        rule.mapTo(domIt.lowerCorner(), domIt.upperCorner(),nodes,weights);
        damage.piece(0).eval_into(nodes,vals);
        labels[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
    }

    return labels;
}

template <short_t dim, class T>
T refineMesh(         gsMultiBasis<T>& basis,
                const std::vector<T> & vals,
                const gsOptionList   & options)
{
    typedef typename gsHElementHelper<dim,T>::HElementContainer HElementContainer;

    gsHElementMarker<dim,T> marker(basis.basis(0));
    marker.options().setSwitch("Admissible",true);
    marker.options().setInt("MaxLevel",1);
    marker.options().setInt("RefineRule",1);
    marker.options().setInt("CoarsenRule",1);
    marker.options().setReal("RefineParam",0.1);
    marker.options().setReal("CoarsenParam",0.1);
    marker.options().update(options,gsOptionList::ignoreIfUnknown);

    marker.setErrors(vals);
    HElementContainer markedRef = marker.markRef();

    T area = 0.0;
    gsMatrix<T> box;
    for (const auto & elem : markedRef)
    {
        box = marker.helper().toBox(elem);
        area += (box.col(1)-box.col(0)).prod();
    }

    std::vector<index_t> refBox = marker.toRefBoxes(markedRef);
    basis.basis(0).refineElements(refBox);
    return area;
}

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & solverParameters,
           gsOptionList & mesherOptions,
           gsMultiPatch<T> & mp_ini,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           std::string & outputdir)
{
    ////////////////////////////////////////////////////////////////////////////////////
    // Load parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Young's modulus [N/mm^2]
    T E = materialParameters.getReal("E");
    // Poisson's ratio [-]
    T nu = materialParameters.getReal("nu");
    // Toughness [N/mm]
    T Gc = materialParameters.getReal("Gc");
    // Internal length [mm]
    T l0 = materialParameters.getReal("l0");
    // Order of the phase-field model
    index_t order = materialParameters.getInt("order");
    // AT model
    index_t AT = materialParameters.getInt("AT");

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    ////////////////////////////////////////////////////////////////////////////////////
    // Boundary control parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Maximum displacement [mm]
    T uend = controlParameters.getReal("uend");
    // Displacement step [mm]
    T ustep = controlParameters.getReal("ustep");
    // Initial displacement [mm]
    T ucurr = controlParameters.getReal("umin");
    // Step transition [mm]
    T utrans = controlParameters.askReal("utrans",uend);
    // Step reduction factor [-]
    T ured = controlParameters.askReal("ured",1.);
    // Maximum number of iterations
    index_t maxIt = controlParameters.getInt("maxIt");
    // Maximum number of iterations for elasticity problem
    index_t maxItEl = controlParameters.getInt("maxItEl");
    // Maximum number of iterations for phase-field problem
    index_t maxItPf = controlParameters.getInt("maxItPf");
    // Tolerance for the elasticity problem
    T tolEl = controlParameters.getReal("tolEl");
    // Tolerance for the phase-field problem
    T tolPf = controlParameters.getReal("tolPf");
    // Staggered tolerance
    T tol = controlParameters.getReal("tol");

    // ---- Three-point-bending geometry ratios (relative to the beam length L) ----
    // Left support:  x = supportLeftRatio  * L   (pin:    u_x=u_y=0)
    // Right support: x = supportRightRatio * L   (roller: u_y=0)
    // Loading point: x = loadRatio         * L   (top edge, prescribed u_y = -ucurr)
    T supportLeftRatio  = controlParameters.askReal("supportLeftRatio", 0.1);
    T supportRightRatio = controlParameters.askReal("supportRightRatio",0.9);
    T loadRatio         = controlParameters.askReal("loadRatio",        0.5);
    // In 3D: apply the load as a line (default, spanning the full depth,
    // like a loading roller) or as a single point indenter at z=zMid.
    bool loadAsPoint = controlParameters.askSwitch("loadAsPoint", false);

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Convert to THB
    gsMultiPatch<T> mp;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        // Check if tensor basis
        if      ((dynamic_cast<const gsTensorBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            // Create a THB spline basis
            const gsTensorBSpline<dim,T> & tb = static_cast<const gsTensorBSpline<dim,T> &>(mp_ini.patch(i));
            gsTHBSpline<dim,T> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if ((dynamic_cast<const gsTHBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            const gsTHBSpline<dim,T> & thb = static_cast<const gsTHBSpline<dim,T> &>(mp_ini.patch(i));
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else
            GISMO_ERROR("The basis is not a TB-spline basis or THB-spline basis.");
    }

    // Construct the basis
    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    // Boundary conditions: bc_u is used as-is (typically empty for 3PB); the
    // support/load conditions are imposed separately, see below.
    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);

    // Physical extent of the beam (assumes a straight, axis-aligned rectangle/box)
    gsVector<T> bbLower, bbUpper;
    getBoundingBox<T>(mp,bbLower,bbUpper);
    const T L = bbUpper(0)-bbLower(0);
    const T xLeft  = bbLower(0) + supportLeftRatio *L;
    const T xRight = bbLower(0) + supportRightRatio*L;
    const T xLoad  = bbLower(0) + loadRatio        *L;
    // Center of the depth (z) direction, only meaningful in 3D; used to pick
    // a single "pin" point on the support/load lines.
    const T zMid = (dim==3) ? (T)0.5*(bbLower(2)+bbUpper(2)) : (T)0.0;
    gsInfo<<"Beam length L = "<<L<<"\n";
    gsInfo<<"Left support at  x = "<<xLeft <<"\n";
    gsInfo<<"Right support at x = "<<xRight<<"\n";
    gsInfo<<"Loading point at x = "<<xLoad <<"\n";

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<T> mp_def = mp;
    gsMultiPatch<T> displacement = mp;
    gsMultiPatch<T> ddisplacement;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,damage,dim);
    // Initialize the elasticity assembler
    gsVector<T> bodyForceVec(dim);
    bodyForceVec.setZero();
    gsConstantFunction<T> bodyForce(bodyForceVec,dim);

    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<T> * pfAssembler;

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<T> u, du;

#ifdef GISMO_WITH_PARDISO
    gsInfo<<"Using Pardiso direct solver\n";
#else
    gsInfo<<"Using CG diagonal preconditioned iterative solver\n";
#endif

    times<T> stagTimes;
    times<T> stepTimes;
    stagTimes.reset();
    stepTimes.reset();

    gsSparseMatrix<T> elMatrix;
    gsMatrix<T> elRhs;

    gsMatrix<T> R;
    gsMatrix<T> D, deltaD;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");
    gsParaviewCollection meshCollection(outputdir+"mesh");
    gsStopwatch smallClock, bigClock;

    std::ofstream file;
    file.open(outputdir+"results.txt");
    file<<"LoadStep,u,Fx,Fy,FSupportLeft,FSupportRight,E_u,E_d,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,projectionTime,basis_size,ref_area,totIt_el,totIt_pf,numIt_stag,numIt_ref\n";
    file.close();
    file.open(outputdir+"iteration_results.txt");
    file<<"LoadStep,RefIt,StagIt,u,Unorm,Dnorm,Rnorm,Fnorm,relRnorm,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,basis_size,ref_area,numIt_el,numIt_pf\n";
    file.close();

    T Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    index_t numIt_el = 0, totIt_el = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_pf = 0, totIt_pf = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_ref = 0; // number of refinement iterations per load step
    index_t numIt_stag = 0; // number of staggered iterations per load step

    index_t step = 0;
    while (ucurr<=uend)
    {
        numIt_ref = numIt_stag = 0;
        totIt_el = totIt_pf = 0;
        stepTimes.reset();
        stagTimes.reset();

        bool refined = true;
        index_t basis_size_old, basis_size;
        T basis_size_ratio;
        T markedArea = 0., tmpArea = 0.;
        gsInfo<<"===========================================================================================================================\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n";
        // Refinement iterations
        index_t refIt = 0;
        while(true)
        {
            basis_size = basis_size_old = mb.basis(0).size();
            gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
            gsInfo<<"Refinement iteration "<<refIt<<" (basis size: "<<basis_size<<"):\n";
            // CONSTRUCT ASSEMBLERS (since refreshing is not possible)
            gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
            elAssembler.options().setReal("ExprAssembler.quA",1.0);
            elAssembler.options().setInt ("ExprAssembler.quB",0);
            elAssembler.options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            elAssembler.initialize();
            // Initialize u
            elAssembler.constructSolution(displacement,u);

            // ---- 3PB: locate the support/load dofs for the CURRENT basis ----
            // Locating must be redone every time the basis changes (adaptive
            // refinement), since local dof indices are only meaningful for a
            // fixed basis.
            gsDofMapper mapper = buildCompanionMapper(mb,bc_u,dim);
            std::vector<gsPointConstraint<T>> pdofs = findSupportDofs(mb,mp,xLeft,xRight,zMid,dim);
            std::vector<gsPointConstraint<T>> loadDofs = findLoadDofs(mb,mp,xLoad,ucurr,zMid,dim,loadAsPoint);
            pdofs.insert(pdofs.end(),loadDofs.begin(),loadDofs.end());

            // Helper: (re-)assemble the elasticity operator at u and impose
            // the point constraints on the resulting linear system.
            auto assembleConstrained = [&](const gsMatrix<T> & uvec)
            {
                elAssembler.assemble(uvec);
                elAssembler.matrix_into(elMatrix);
                elAssembler.rhs_into(elRhs);
                applyPointDirichlet(elMatrix,elRhs,mapper,pdofs);
            };

            if      (order == 2 && AT == 1)
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
            else if (order == 4 && AT == 1)
            {
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
                pfAssembler->options().setReal("cw",4.44847);
            }
            else if (order == 2 && AT == 2)
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
            else if (order == 4 && AT == 2)
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
            else
                GISMO_ERROR("Invalid order and/or AT model");

            // Pre-assemble the phase-field operators that do not depend on the solution
            pfAssembler->options().setReal("l0",l0);
            pfAssembler->options().setReal("Gc",Gc);
            pfAssembler->options().setReal("ExprAssembler.quA",1.0);
            pfAssembler->options().setInt ("ExprAssembler.quB",0);
            pfAssembler->options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            pfAssembler->setSpaceBasis(mb);
            pfAssembler->initialize();

            // Construct Phase-Field solution vector from projected solution
            pfAssembler->constructSolution(damage,D);

            // Pre-assemble the phase-field operators that do not depend on the solution
            smallClock.restart();
            pfAssembler->assemblePhi();
            pfAssembler->matrix_into(QPhi);
            pfAssembler->rhs_into(q);
            T pfAssemblyTime0 = smallClock.stop();

            deltaD.setZero(D.rows(),1);
            index_t stagIt = 0;

            while(true)
            {
                stagTimes.reset();
                if (stagIt==0) stagTimes.pfAssemblyTime = pfAssemblyTime0; // save the time of the first assembly

                bigClock.restart();
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<stagIt<<"--------------------------\n";
                gsInfo<<"    ---------------------------------ELASTICITY----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||F||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<PRINT(20)<<"it. solver [s]"<<PRINT(20)<<"solver it."<<"|\n";

                material.setParameter(2,damage);
                // Pre-assemble the elasticity problem
                smallClock.restart();
                assembleConstrained(u);
                stagTimes.elAssemblyTime += smallClock.stop();
                Fnorm = elRhs.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                index_t elIt = 0;
                while(true)
                {
                    // Solve
                    T itSolverTime = 0;
                    index_t itSolverIterations = 0;
                    smallClock.restart();
#ifdef GISMO_WITH_PARDISO
                    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
                    typename gsSparseSolver<T>::CGDiagonal solver;
#endif
                    solver.compute(elMatrix);
                    u = solver.solve(elRhs);
#ifdef GISMO_WITH_PARDISO
                    itSolverIterations = 1;
#else
                    itSolverIterations = solver.iterations();
#endif
                    itSolverTime = smallClock.stop();
                    stagTimes.elSolverTime += itSolverTime;

                    // Re-assemble
                    smallClock.restart();
                    assembleConstrained(u);
                    stagTimes.elAssemblyTime += smallClock.stop();
                    Fnorm = elRhs.norm();
                    Fnorm = (Fnorm == 0) ? 1 : Fnorm;

                    // Check convergence with the old matrix and rhs (saves one assembly)
                    Rnorm = (elMatrix*u - elRhs).norm();
                    gsInfo<<"    | "<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm<<PRINT(14)<<Fnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<stagTimes.elAssemblyTime<<PRINT(20)<<stagTimes.elSolverTime<<PRINT(20)<<itSolverTime<<PRINT(20)<<itSolverIterations<<"|\n";

                    if (Rnorm/Fnorm < tolEl || u.norm() < 1e-12 || maxItEl==1)
                        break;
                    else if (elIt == maxItEl-1)
                        GISMO_ERROR("Elasticity problem did not converge.");
                    else
                        elIt++;
                }

                numIt_el = elIt+1;
                totIt_el+= numIt_el;

                elAssembler.constructSolution(u,displacement);
                for (size_t p=0; p!=mp.nPatches(); ++p)
                    mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

                // Initialize the function for the elastic energy
                gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

                // ==================================================================================

                gsInfo<<"    ---------------------------------PHASE-FIELD---------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"|\n";

                // Phase-field problem
                smallClock.restart();
                pfAssembler->assemblePsi(Psi);
                stagTimes.pfAssemblyTime += smallClock.stop();
                pfAssembler->matrix_into(QPsi);
                pfAssembler->rhs_into(qpsi);
                if (qpsi.rows()==0) // qpsi is empty for AT2 models
                    qpsi = gsMatrix<T>::Zero(QPsi.rows(),1);
                Q = QPhi + QPsi;

                // Reconstruct the solution from the damage field
                pfAssembler->constructSolution(damage,D);

                // Initialize the PSOR solver
                smallClock.restart();
                gsPSOR<T> PSORsolver(Q);
                PSORsolver.options().setInt("MaxIterations",30000);
                PSORsolver.options().setSwitch("Verbose",false);
                PSORsolver.options().setReal("tolU",1e-4);
                PSORsolver.options().setReal("tolNeg",1e-9);
                PSORsolver.options().setReal("tolPos",1e-9);
                stagTimes.pfSolverTime = smallClock.stop();
                index_t pfIt = 0;
                while(true)
                {
                    // Assemble
                    smallClock.restart();
                    R = Q * D - qpsi + q;
                    stagTimes.pfAssemblyTime += smallClock.stop();

                    // Solve
                    smallClock.restart();
                    PSORsolver.solve(R,deltaD); // deltaD = Q \ R
                    stagTimes.pfSolverTime += smallClock.stop();
                    D += deltaD;
                    gsInfo<<"    | "<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<stagTimes.pfAssemblyTime<<PRINT(20)<<stagTimes.pfSolverTime<<"|\n";
                    if (deltaD.norm()/D.norm() < tolPf || D.norm() < 1e-12 || maxItPf==1)
                        break;
                    else if (pfIt == maxItPf-1 && maxItPf != 1)
                        GISMO_ERROR("Phase-field problem did not converge.");
                    pfIt++;
                }
                numIt_pf = pfIt+1;
                totIt_pf+= numIt_pf;

                // Update damage spline
                pfAssembler->constructSolution(D,damage);

                material.setParameter(2,damage);
                smallClock.restart();
                assembleConstrained(u);
                stagTimes.elAssemblyTime += smallClock.stop();
                Fnorm = elRhs.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                Rnorm = (elMatrix*u - elRhs).norm();

                gsInfo<<"    ----------------------------------FINISHED-----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"||R|| = "<<PRINT(14)<<Rnorm<<PRINT(20)<<"total time [s] = "<<PRINT(20)<<bigClock.stop()<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||F|| = "<<PRINT(14)<<Fnorm<<PRINT(20)<<"elasticity time [s] = "<<PRINT(20)<<stagTimes.elAssemblyTime+stagTimes.elSolverTime<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||D|| = "<<PRINT(14)<<D.norm()<<PRINT(20)<<"phase-field time [s] = "<<PRINT(20)<<stagTimes.pfAssemblyTime+stagTimes.pfSolverTime<<"\n";
                gsInfo<<"    -----------------------------------------------------------------------------\n";

                file.open(outputdir+"iteration_results.txt",std::ios::app);
                file<<step<<","<<refIt<<","<<stagIt<<","<<ucurr<<","
                    <<u.norm()<<","<<D.norm()<<","<<Rnorm<<","<<Fnorm<<","<<Rnorm/Fnorm<<","
                    <<stagTimes.elAssemblyTime<<","<<stagTimes.elSolverTime<<","
                    <<stagTimes.pfAssemblyTime<<","<<stagTimes.pfSolverTime<<","
                    <<basis_size<<","<<markedArea<<","
                    <<numIt_el<<","<<numIt_pf<<"\n";
                file.close();

                stepTimes.elAssemblyTime += stagTimes.elAssemblyTime;
                stepTimes.elSolverTime += stagTimes.elSolverTime;
                stepTimes.pfAssemblyTime += stagTimes.pfAssemblyTime;
                stepTimes.pfSolverTime += stagTimes.pfSolverTime;

                if (Rnorm/Fnorm < tol)
                    break;
                else if (stagIt == maxIt-1)
                    GISMO_ERROR("Staggered iterations problem did not converge.");
                else
                    stagIt++;
            }
            numIt_stag += stagIt+1;

            gsInfo<<"\n";
            gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";

            // =========================================================================
            // REFINE MESH
            // All labelled elements are refined to the maximum level, step-by-step
            for (index_t i=0; i!=mesherOptions.askInt("MaxLevel",1); ++i)
            {
                smallClock.restart();
                std::vector<T> elVals = labelElements<dim,T>(mp, damage, mb,0.1,1.0);
                gsInfo<<"Labelling level "<<i<<" took "<<smallClock.stop()<<" seconds\n";
                if (gsAsVector<T>(elVals).sum() > 0)
                {
                    smallClock.restart();
                    tmpArea = refineMesh<dim,T>(mb,elVals,mesherOptions);
                    gsInfo<<"Refining mesh took "<<smallClock.stop()<<" seconds\n";
                }

                tmpArea /= (mb.basis(0).support().col(1)-mb.basis(0).support().col(0)).prod();
                markedArea = math::max(markedArea,tmpArea);
                basis_size = mb.basis(0).size();
                refined = basis_size > basis_size_old;
                if (!refined)
                    break;
            }
            gsInfo<<"Marked area: "<<markedArea<<"\n";
            refined &= markedArea > mesherOptions.askReal("SizeRatio",1.01);

            // =========================================================================
            // PROJECT SOLUTIONS
            smallClock.restart();
            gsMatrix<T> projCoefs;
            // Geometry
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
            mp.clear();
            mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            mp_def = mp;
            // Displacement
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
            displacement.clear();
            displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // Damage
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),damage.patch(0),projCoefs);
            damage.clear();
            damage.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            T ptime = smallClock.stop();
            stepTimes.projectionTime += ptime;
            gsInfo<<"Projection took "<<ptime<<" seconds\n";

            delete pfAssembler;

            if (!refined)
                break;
            else if (refIt==mesherOptions.askInt("MaxRefIterations",5)-1)
            {
                gsWarn<<"Maximum number of refinement itertions reached\n";
                break;
            }
            else
                refIt++;
        }
        numIt_ref = refIt+1;

        // =========================================================================
        // PLOT
        if (plot && step%plotmod==0)
        {
            std::string filename;

            filename = "mesh_" + util::to_string(step);
            gsMesh<T> mesh(mb.basis(0));
            writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename,1);
            meshCollection.addPart(filename,step,"Mesh",0);

            filename = "damage_" + util::to_string(step);
            gsWriteParaview(mp,damage,outputdir+filename,100000);
            filename += "0";
            damageCollection.addPart(filename,step,"Solution",0);

            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
            filename = "Psi_"+util::to_string(step);
            gsWriteParaview(mp,Psi,outputdir+filename,100000);
            filename += "0";
            psiCollection.addPart(filename,step,"Solution",0);

            filename = "displacement_"+util::to_string(step);
            gsWriteParaview(mp,displacement,outputdir+filename,1000);
            filename += "0";
            displCollection.addPart(filename,step,"Solution",0);
        }

        // =========================================================================
        // REACTION FORCES at the load point and at the two supports.
        //
        // fullElAssembler uses an EMPTY gsBoundaryConditions object, so ALL dofs
        // are free and its internal numbering coincides with a plain,
        // non-eliminating gsDofMapper(mb,dim) -- exactly the trick already used
        // in the original example to compute side reactions, here reused to
        // read out reactions at individual (support/load) dofs instead of a
        // whole side.
        gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
        fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
        fullElAssembler.options().setInt ("ExprAssembler.quB",0);
        fullElAssembler.initialize();

        // Compute resulting force and energies
        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull);
        gsMatrix<T> Rfull = fullElAssembler.matrix()*ufull - fullElAssembler.rhs();
        gsDofMapper fullMapper(mb,dim);
        fullMapper.finalize();

        // Sum the reaction over the whole support/load line (in 2D each line
        // has a single dof, so this reduces to the original single-dof read-out).
        gsVector<T> loc;
        const std::vector<index_t> loadLine         = locateBoundaryLineAtX(mb,mp,0,boundary::north,xLoad, loc);
        const std::vector<index_t> supportLeftLine  = locateBoundaryLineAtX(mb,mp,0,boundary::south,xLeft, loc);
        const std::vector<index_t> supportRightLine = locateBoundaryLineAtX(mb,mp,0,boundary::south,xRight,loc);

        T Fx = 0, Fy = 0, FSupportLeft = 0, FSupportRight = 0;
        for (const index_t dof : loadLine)
        {
            Fx += Rfull(fullMapper.index(dof,0,0),0);
            Fy += Rfull(fullMapper.index(dof,0,1),0);
        }
        for (const index_t dof : supportLeftLine)
            FSupportLeft  += Rfull(fullMapper.index(dof,0,1),0);
        for (const index_t dof : supportRightLine)
            FSupportRight += Rfull(fullMapper.index(dof,0,1),0);

        file.open(outputdir+"results.txt",std::ios::app);
        file<<step<<","<<ucurr<<","<<-Fx<<","<<-Fy<<","<<-FSupportLeft<<","<<-FSupportRight<<","
            <<(0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value()<<","
            <<(0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value()<<","
            <<stepTimes.elAssemblyTime<<","<<stepTimes.elSolverTime<<","
            <<stepTimes.pfAssemblyTime<<","<<stepTimes.pfSolverTime<<","
            <<stepTimes.projectionTime<<","
            <<basis_size<<","<<markedArea<<","
            <<totIt_el<<","<<totIt_pf<<","
            <<numIt_stag<<","<<numIt_ref<<"\n";
        file.close();

        // =========================================================================
        // INCREMENT STEP

        if (ucurr == uend || math::abs(ucurr-uend) < 1e-10)
            break;

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;
        ucurr = math::min(ucurr,uend);
        step++;
    }

    if (plot)
    {
        meshCollection.save();
        damageCollection.save();
        psiCollection.save();
        displCollection.save();
    }
}
