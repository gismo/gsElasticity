/** @file adaptive_phase-field_fracture_example_DC_3pb_parametric.cpp

    @brief Three-point bending phase-field fracture example with point
    constraints applied in PARAMETRIC coordinates.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    -----------------------------------------------------------------------
    Variant of adaptive_phase-field_fracture_example_DC_3pb.cpp.

    The three-point-bending support/load constraints are imposed at
    PARAMETRIC (not physical) coordinates using active_into + eval_into
    to form multi-DOF constraints, then applying them by constrained
    elimination (applyConstraint). This replaces the physical-coordinate
    DOF search (locateBoundaryLineAtX / findSupportDofs / findLoadDofs)
    and the single-DOF applyPointDirichlet with a single helper that
    handles higher-order, non-interpolatory bases correctly.

    Boundary conditions (xy is the bending plane, z is out-of-plane):
      Left  roller (parametric x = supportLeftRatio,  y = 0): u_y = 0, u_z = 0 (3D)
      Right roller (parametric x = supportRightRatio, y = 0): u_y = 0, u_z = 0 (3D)
      Load         (parametric x = loadRatio,         y = 1): u_x = 0, u_z = 0 (3D),
                                                               u_y = ucurr (prescribed)

    In 3D each line constraint is enforced by calling applyConstraint once
    per unique z-Greville abscissa found on the relevant boundary side of
    the current (possibly THB-refined) basis.

    To run in 3D:
    ./bin/adaptive_phase-field_fracture_example_DC_3pb_parametric \
        -I optional/gsElasticity/filedata/phase-field-fracture/AT-1_Order4/3PB_THB_3D --plot
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
// 3-POINT-BENDING PARAMETRIC CONSTRAINT HELPERS //////////////////////////
//////////////////////////////////////////////////////////////////////////

/** \brief Impose a scalar constraint c^T u_{comp} = value at a parametric
 *  point \a pt by exact algebraic constrained elimination on the assembled
 *  linear system (K, F).
 *
 *  The constraint vector c_i = N_i(pt) is built from the basis functions
 *  active at \a pt via active_into / eval_into.  One "master" degree of
 *  freedom (that with the largest |c_i|) is eliminated in favour of the
 *  others, following the same Dirichlet-elimination pattern as the
 *  single-DOF applyPointDirichlet but extended to multi-DOF constraints:
 *
 *    F[r]         -= K[r,m] * value / c_m     for r != m   (same as single-DOF)
 *    K[r, gd[i]]  -= K[r,m] * c[i]  / c_m    for r != m, i != master  (new term)
 *
 *  After the substitution the constraint row is set to c^T u_{comp} = value,
 *  so the solved vector directly contains the constrained displacement.
 *
 *  For a clamped B-spline at a boundary point (y=0 or y=1), all active
 *  basis functions in the transverse direction collapse to the first (or
 *  last) row, so the constraint only involves the in-plane DOFs at that
 *  boundary -- exactly the desired behaviour for a cylindrical roller/load.
 */
template <class T>
void applyConstraint(gsSparseMatrix<T> & K,
                      gsMatrix<T>       & F,
                      const gsMultiBasis<T> & mb,
                      const gsDofMapper & mapper,
                      const gsMatrix<T> & pt,     // parametric point (dim x 1)
                      index_t             comp,   // 0=u_x, 1=u_y, 2=u_z
                      T                   value)  // prescribed value
{
    gsMatrix<index_t> active;
    gsMatrix<T>       vals;
    mb.basis(0).active_into(pt, active);   // active: (na x 1)
    mb.basis(0).eval_into  (pt, vals);     // vals:   (1  x na)

    const index_t na = active.rows();
    std::vector<index_t> gd(na);
    std::vector<T>       c(na);
    for (index_t i = 0; i < na; ++i)
    {
        gd[i] = mapper.index(active(i,0), 0, comp);
        c[i]  = vals(0,i);
    }

    // Master DOF: largest |c_i| -- best conditioned for elimination
    index_t best = 0;
    for (index_t i = 1; i < na; ++i)
        if (math::abs(c[i]) > math::abs(c[best])) best = i;
    const index_t m  = gd[best];
    const T       cm = c[best];

    // Step 1: column-m substitution into all rows r != m.
    // Collect column-m entries first to avoid iterator invalidation when
    // coeffRef inserts new entries into the sparse matrix.
    std::vector<std::pair<index_t,T>> col_m;
    for (typename gsSparseMatrix<T>::InnerIterator it(K, m); it; ++it)
        if (it.index() != m)
            col_m.emplace_back(it.index(), it.value());

    for (const auto & entry : col_m)
    {
        const index_t r    = entry.first;
        const T       k_rm = entry.second;
        F(r, 0) -= k_rm * value / cm;
        for (index_t i = 0; i < na; ++i)
            if (gd[i] != m)
                K.coeffRef(r, gd[i]) -= k_rm * c[i] / cm;
    }

    // Step 2: zero row m and column m (same pattern as applyPointDirichlet)
    for (index_t k = 0; k < K.outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(K, k); it; ++it)
            if (it.index() == m || k == m)
                it.valueRef() = T(0);

    // Step 3: replace row m with the constraint equation c^T u_{comp} = value
    for (index_t i = 0; i < na; ++i)
        K.coeffRef(m, gd[i]) = c[i];
    F(m, 0) = value;
}

/** \brief Build a gsDofMapper whose numbering matches the internal numbering
 *  used by a gsSolidAssembler constructed with \a bc.
 *  (Copied verbatim from adaptive_phase-field_fracture_example_DC_3pb.cpp.)
 */
template <class T>
gsDofMapper buildCompanionMapper(const gsMultiBasis<T>         & mb,
                                  const gsBoundaryConditions<T> & bc,
                                  index_t                         dim)
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

/** \brief Return the list of parametric points at which to enforce the
 *  line constraint at (xParam, yParam) in the z-direction.
 *
 *  In 2D: returns the single point {xParam, yParam}.
 *  In 3D: returns one point per unique z-Greville abscissa found among
 *  the boundary DOFs on \a bndSide, so that the full z-line is exactly
 *  constrained (one independent elimination per z-Greville point).
 */
template <short_t dim, class T>
std::vector<gsMatrix<T>> constraintPoints(
    const gsMultiBasis<T> & mb,
    boxSide                 bndSide,
    T                       xParam,
    T                       yParam)
{
    std::vector<gsMatrix<T>> pts;
    if (dim == 2)
    {
        gsMatrix<T> pt(2,1);
        pt << xParam, yParam;
        pts.push_back(pt);
    }
    else // dim == 3
    {
        gsMatrix<index_t> bnd     = mb.basis(0).boundary(bndSide);
        gsMatrix<T>       anchors = mb.basis(0).anchors();   // dim x n_dofs

        // Collect unique z-anchor values (parametric component dim-1)
        std::vector<T> z_vals;
        for (index_t k = 0; k < bnd.rows(); ++k)
        {
            const T z = anchors(2, bnd(k,0));
            bool found = false;
            for (const T & zv : z_vals)
                if (math::abs(z - zv) < (T)1e-10) { found = true; break; }
            if (!found) z_vals.push_back(z);
        }
        std::sort(z_vals.begin(), z_vals.end());

        for (const T z_k : z_vals)
        {
            gsMatrix<T> pt(3,1);
            pt << xParam, yParam, z_k;
            pts.push_back(pt);
        }
    }
    return pts;
}

/** \brief Find all local dofs on a boundary side whose first parametric
 *  component (x-anchor) is closest to \a paramX.
 *
 *  Parametric-coordinate counterpart of locateBoundaryLineAtX: works
 *  directly in parameter space without geometry evaluation.
 */
template <class T>
std::vector<index_t> findBoundaryLineByParamX(
    const gsMultiBasis<T> & mb,
    boxSide                 side,
    T                       paramX,
    T                       tol = (T)1e-6)
{
    gsMatrix<index_t> bnd     = mb.basis(0).boundary(side);
    gsMatrix<T>       anchors = mb.basis(0).anchors();

    T bestDist = std::numeric_limits<T>::max();
    for (index_t k = 0; k < bnd.rows(); ++k)
    {
        const T d = math::abs(anchors(0, bnd(k,0)) - paramX);
        if (d < bestDist) bestDist = d;
    }

    std::vector<index_t> result;
    for (index_t k = 0; k < bnd.rows(); ++k)
        if (math::abs(anchors(0, bnd(k,0)) - paramX) <= bestDist + tol)
            result.push_back(bnd(k,0));
    return result;
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

    gsCmdLine cmd("Three-point bending phase-field fracture (parametric point constraints).");
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
    std::string parInputPath    = (parInput.empty()    ? inputDir + "parameters.xml" : parInput);
    std::string geoInputPath    = (geoInput.empty()    ? inputDir + "geometry.xml"   : geoInput);
    std::string damageInputPath = (damageInput.empty() ? inputDir + "damage.xml"     : damageInput);
    GISMO_ASSERT(gsFileManager::fileExists(parInputPath),    "Input parameter file "<<parInputPath   <<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(geoInputPath),    "Input geometry file "  <<geoInputPath   <<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(damageInputPath), "Input damage file "    <<damageInputPath<<" not found.");
    gsInfo << "Input parameter file "<< parInputPath    <<"\n";
    gsInfo << "Input geometry file " << geoInputPath    <<"\n";
    gsInfo << "Input damage file "   << damageInputPath <<"\n";

    if (output.empty()) output = "./output/";
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
    GISMO_ASSERT(fd_pars.hasLabel("control"),  "Displacement-Control parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("meshing"),  "Adaptive meshing parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_u"),    "Displacement boundary conditions not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_d"),    "Phase-field boundary conditions not found in the input file.");

    gsOptionList materialParameters;
    fd_pars.getLabel("material", materialParameters);

    gsOptionList controlParameters;
    fd_pars.getLabel("control", controlParameters);

    gsOptionList solverParameters;
    if (fd_pars.hasLabel("solver"))
        fd_pars.getLabel("solver", solverParameters);

    // NOTE: for a pure 3-point-bending setup, BCs_u is typically EMPTY:
    // support/load constraints are applied as parametric point constraints.
    // Additional side conditions (e.g. out-of-plane symmetry) remain possible.
    gsBoundaryConditions<> bc_u;
    fd_pars.getLabel("BCs_u", bc_u);

    gsBoundaryConditions<> bc_d;
    fd_pars.getLabel("BCs_d", bc_d);

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);

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
    gsLobattoRule<T> rule(np);
    for (typename gsBasis<T>::domainIter domIt = basis.basis(0).domain()->beginAll(); domIt<domEnd; ++domIt)
    {
        gsMatrix<T> nodes, vals;
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
    T E      = materialParameters.getReal("E");
    T nu     = materialParameters.getReal("nu");
    T Gc     = materialParameters.getReal("Gc");
    T l0     = materialParameters.getReal("l0");
    index_t order = materialParameters.getInt("order");
    index_t AT    = materialParameters.getInt("AT");

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2,       "Please specify the AT model (1 or 2).");

    ////////////////////////////////////////////////////////////////////////////////////
    // Boundary control parameters
    ////////////////////////////////////////////////////////////////////////////////////
    T uend   = controlParameters.getReal("uend");
    T ustep  = controlParameters.getReal("ustep");
    T ucurr  = controlParameters.getReal("umin");
    T utrans = controlParameters.askReal("utrans",uend);
    T ured   = controlParameters.askReal("ured",1.);
    index_t maxIt   = controlParameters.getInt("maxIt");
    index_t maxItEl = controlParameters.getInt("maxItEl");
    index_t maxItPf = controlParameters.getInt("maxItPf");
    T tolEl  = controlParameters.getReal("tolEl");
    T tolPf  = controlParameters.getReal("tolPf");
    T tol    = controlParameters.getReal("tol");

    // Parametric ratios for support / load positions (domain [0,1])
    T supportLeftRatio  = controlParameters.askReal("supportLeftRatio", 0.1);
    T supportRightRatio = controlParameters.askReal("supportRightRatio",0.9);
    T loadRatio         = controlParameters.askReal("loadRatio",        0.5);
    // In 3D: apply load as a line spanning full z (default) or as a point at z=0.5
    bool loadAsPoint = controlParameters.askSwitch("loadAsPoint", false);

    gsInfo<<"Left support at  parametric x = "<<supportLeftRatio <<"\n";
    gsInfo<<"Right support at parametric x = "<<supportRightRatio<<"\n";
    gsInfo<<"Loading point at parametric x = "<<loadRatio        <<"\n";

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Convert to THB (required for adaptive refinement)
    gsMultiPatch<T> mp;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        if (dynamic_cast<const gsTensorBSpline<dim,T> *>(&mp_ini.patch(i)))
        {
            const gsTensorBSpline<dim,T> & tb = static_cast<const gsTensorBSpline<dim,T> &>(mp_ini.patch(i));
            gsTHBSpline<dim,T> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if (dynamic_cast<const gsTHBSpline<dim,T> *>(&mp_ini.patch(i)))
        {
            const gsTHBSpline<dim,T> & thb = static_cast<const gsTHBSpline<dim,T> &>(mp_ini.patch(i));
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else
            GISMO_ERROR("The basis is not a TensorBSpline or THBSpline.");
    }

    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    gsMultiPatch<T> mp_def = mp;
    gsMultiPatch<T> displacement = mp;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    gsLinearDegradedMaterial<T> material(E,nu,damage,dim);

    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    gsPhaseFieldAssemblerBase<T> * pfAssembler;

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<T> u;

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
    index_t numIt_el = 0, totIt_el = 0;
    index_t numIt_pf = 0, totIt_pf = 0;
    index_t numIt_ref = 0;
    index_t numIt_stag = 0;

    index_t step = 0;
    while (ucurr<=uend)
    {
        numIt_ref = numIt_stag = 0;
        totIt_el = totIt_pf = 0;
        stepTimes.reset();
        stagTimes.reset();

        bool refined = true;
        index_t basis_size_old, basis_size;
        T markedArea = 0., tmpArea = 0.;
        gsInfo<<"===========================================================================================================================\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n";
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
            elAssembler.constructSolution(displacement,u);

            // ---- 3PB: build parametric constraint points for the CURRENT basis ----
            // Must be redone every refinement iteration since z-Greville points
            // may change as the THB basis is hierarchically refined.
            gsDofMapper mapper = buildCompanionMapper(mb,bc_u,dim);

            std::vector<gsMatrix<T>> pts_left  = constraintPoints<dim,T>(mb, boundary::south, supportLeftRatio,  (T)0);
            std::vector<gsMatrix<T>> pts_right = constraintPoints<dim,T>(mb, boundary::south, supportRightRatio, (T)0);
            std::vector<gsMatrix<T>> pts_load;
            if (dim == 2)
            {
                gsMatrix<T> pt(2,1); pt << loadRatio, (T)1;
                pts_load.push_back(pt);
            }
            else // dim == 3
            {
                if (loadAsPoint)
                {
                    gsMatrix<T> pt(3,1); pt << loadRatio, (T)1, (T)0.5;
                    pts_load.push_back(pt);
                }
                else
                    pts_load = constraintPoints<dim,T>(mb, boundary::north, loadRatio, (T)1);
            }

            // Helper: (re-)assemble K,F and apply all 3PB parametric constraints.
            // Rollers (y=0): u_y=0  and (3D) u_z=0
            // Load    (y=1): u_x=0, u_y=ucurr, and (3D) u_z=0
            auto assembleConstrained = [&](const gsMatrix<T> & uvec)
            {
                elAssembler.assemble(uvec);
                elAssembler.matrix_into(elMatrix);
                elAssembler.rhs_into(elRhs);

                for (const auto & pt : pts_left)
                {
                    applyConstraint(elMatrix, elRhs, mb, mapper, pt, 1, (T)0);   // u_y = 0
                    if (dim == 3)
                        applyConstraint(elMatrix, elRhs, mb, mapper, pt, 2, (T)0); // u_z = 0
                }
                for (const auto & pt : pts_right)
                {
                    applyConstraint(elMatrix, elRhs, mb, mapper, pt, 1, (T)0);   // u_y = 0
                    if (dim == 3)
                        applyConstraint(elMatrix, elRhs, mb, mapper, pt, 2, (T)0); // u_z = 0
                }
                for (const auto & pt : pts_load)
                {
                    applyConstraint(elMatrix, elRhs, mb, mapper, pt, 0, (T)0);   // u_x = 0
                    if (dim == 3)
                        applyConstraint(elMatrix, elRhs, mb, mapper, pt, 2, (T)0); // u_z = 0
                    applyConstraint(elMatrix, elRhs, mb, mapper, pt, 1, ucurr);  // u_y = ucurr
                }
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

            pfAssembler->options().setReal("l0",l0);
            pfAssembler->options().setReal("Gc",Gc);
            pfAssembler->options().setReal("ExprAssembler.quA",1.0);
            pfAssembler->options().setInt ("ExprAssembler.quB",0);
            pfAssembler->options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            pfAssembler->setSpaceBasis(mb);
            pfAssembler->initialize();

            pfAssembler->constructSolution(damage,D);

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
                if (stagIt==0) stagTimes.pfAssemblyTime = pfAssemblyTime0;

                bigClock.restart();
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<stagIt<<"--------------------------\n";
                gsInfo<<"    ---------------------------------ELASTICITY----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||F||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<PRINT(20)<<"it. solver [s]"<<PRINT(20)<<"solver it."<<"|\n";

                material.setParameter(2,damage);
                smallClock.restart();
                assembleConstrained(u);
                stagTimes.elAssemblyTime += smallClock.stop();
                Fnorm = elRhs.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                index_t elIt = 0;
                while(true)
                {
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

                    smallClock.restart();
                    assembleConstrained(u);
                    stagTimes.elAssemblyTime += smallClock.stop();
                    Fnorm = elRhs.norm();
                    Fnorm = (Fnorm == 0) ? 1 : Fnorm;

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

                gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

                gsInfo<<"    ---------------------------------PHASE-FIELD---------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"|\n";

                smallClock.restart();
                pfAssembler->assemblePsi(Psi);
                stagTimes.pfAssemblyTime += smallClock.stop();
                pfAssembler->matrix_into(QPsi);
                pfAssembler->rhs_into(qpsi);
                if (qpsi.rows()==0)
                    qpsi = gsMatrix<T>::Zero(QPsi.rows(),1);
                Q = QPhi + QPsi;

                pfAssembler->constructSolution(damage,D);

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
                    smallClock.restart();
                    R = Q * D - qpsi + q;
                    stagTimes.pfAssemblyTime += smallClock.stop();

                    smallClock.restart();
                    PSORsolver.solve(R,deltaD);
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

                pfAssembler->constructSolution(D,damage);

                material.setParameter(2,damage);
                smallClock.restart();
                assembleConstrained(u);
                stagTimes.elAssemblyTime += smallClock.stop();
                Fnorm = elRhs.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                Rnorm = (elMatrix*u - elRhs).norm();

                gsInfo<<"    ----------------------------------FINISHED-----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"||R|| = "<<PRINT(14)<<Rnorm<<PRINT(20)<<"total time [s] = "       <<PRINT(20)<<bigClock.stop()<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||F|| = "<<PRINT(14)<<Fnorm<<PRINT(20)<<"elasticity time [s] = "  <<PRINT(20)<<stagTimes.elAssemblyTime+stagTimes.elSolverTime<<"\n";
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
                stepTimes.elSolverTime   += stagTimes.elSolverTime;
                stepTimes.pfAssemblyTime += stagTimes.pfAssemblyTime;
                stepTimes.pfSolverTime   += stagTimes.pfSolverTime;

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
                if (!refined) break;
            }
            gsInfo<<"Marked area: "<<markedArea<<"\n";
            refined &= markedArea > mesherOptions.askReal("SizeRatio",1.01);

            // =========================================================================
            // PROJECT SOLUTIONS
            smallClock.restart();
            gsMatrix<T> projCoefs;
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
            mp.clear();
            mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            mp_def = mp;
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
            displacement.clear();
            displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
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
                gsWarn<<"Maximum number of refinement iterations reached\n";
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
        // REACTION FORCES
        // fullElAssembler uses empty bc so all DOFs are free and its mapper
        // coincides with gsDofMapper(mb,dim) -- used to read reaction forces
        // at the support/load lines in parametric coordinates.
        gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
        fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
        fullElAssembler.options().setInt ("ExprAssembler.quB",0);
        fullElAssembler.initialize();

        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull);
        gsMatrix<T> Rfull = fullElAssembler.matrix()*ufull - fullElAssembler.rhs();
        gsDofMapper fullMapper(mb,dim);
        fullMapper.finalize();

        // Find boundary DOFs at each constraint line using parametric x
        const std::vector<index_t> loadLine         = findBoundaryLineByParamX(mb, boundary::north, loadRatio);
        const std::vector<index_t> supportLeftLine  = findBoundaryLineByParamX(mb, boundary::south, supportLeftRatio);
        const std::vector<index_t> supportRightLine = findBoundaryLineByParamX(mb, boundary::south, supportRightRatio);

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
