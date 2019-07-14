/** @file gsElasticityAssembler.hpp

    @brief Provides linear and nonlinear elasticity systems for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
#include <gsElasticity/gsVisitorMixedLinearElasticity.h>
#include <gsElasticity/gsVisitorMixedNonLinearElasticity.h>
#include <gsElasticity/gsVisitorNonLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(const gsMultiPatch<T> & patches,
                                                const gsMultiBasis<T> & basis,
                                                const gsBoundaryConditions<T> & bconditions,
                                                const gsFunction<T> & body_force)
{
    // Originally concieved as a meaningful class, now gsPde is just a container for
    // the domain, boundary conditions and the right-hand side;
    // any derived class can surve this purpuse, for example gsPoissonPde;
    // TUDO: change/remove gsPde from gsAssembler logic
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    // gsAssembler<>::initialize requires a vector of bases, one for each unknown;
    // different bases are used to compute Dirichlet DoFs;
    // but always the first basis is used for the assembly;
    // TODO: change gsAssembler logic
    m_dim = body_force.targetDim();
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basis);

    Base::initialize(pde, m_bases, defaultOptions());
}

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
                                                gsMultiBasis<T> const & basisDisp,
                                                gsMultiBasis<T> const & basisPres,
                                                gsBoundaryConditions<T> const & bconditions,
                                                const gsFunction<T> & body_force)
{
    // same as above
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    // same as above
    m_dim = body_force.targetDim();
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basisDisp);
    m_bases.push_back(basisPres);

    Base::initialize(pde, m_bases, defaultOptions());
    m_options.setInt("MaterialLaw",material_law::neo_hooke_ln);
}

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addReal("DirichletAssembly","Dirichlet BC scaling parameter for assembly",1.);
    opt.addReal("DirichletConstruction","Dirichlet BC scaling parameter for solution construction",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

template <class T>
void gsElasticityAssembler<T>::refresh()
{
    GISMO_ENSURE(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ENSURE(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_bases.size());
    for (unsigned d = 0; d < m_bases.size(); d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_bases.size());
    m_system = gsSparseSystem<T>(m_dofMappers, dims);

    m_options.setReal("bdO",m_bases.size()*(1+m_options.getReal("bdO"))-1);
    m_system.reserve(m_bases[0], m_options, 1);

    for (unsigned d = 0; d < m_bases.size(); ++d)
        Base::computeDirichletDofs(d);
}

template<class T>
void gsElasticityAssembler<T>::assemble(bool assembleMatrix)
{
    if (assembleMatrix)
    {
        m_system.matrix().setZero();
        m_system.reserve(m_bases[0], m_options, 1);
    }
    m_system.rhs().setZero(Base::numDofs(),1);

    if ( this->numDofs() == 0 )
    {
        gsWarn << "No internal DOFs. Computed Dirichlet boundary only.\n";
        return;
    }

    scaleDDoFs(m_options.getReal("DirichletAssembly"));

    // Compute volumetric integrals and write to the global linear system
    if (m_bases.size() == unsigned(m_dim)) // displacement formulation
    {
        gsVisitorLinearElasticity<T> visitor(*m_pde_ptr,assembleMatrix);
        Base::template push<gsVisitorLinearElasticity<T> >(visitor);
    }
    else // mixed formulation (displacement + pressure)
    {
        gsVisitorMixedLinearElasticity<T> visitor(*m_pde_ptr,assembleMatrix);
        Base::template push<gsVisitorMixedLinearElasticity<T> >(visitor);
    }

    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
    resetDDoFs();
}

template <class T>
bool gsElasticityAssembler<T>::assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix)
{
    if (m_bases.size() == unsigned(m_dim)) // displacement formulation
    {
        gsMultiPatch<T> displacement;
        constructSolution(solutionVector,displacement);
        if (checkSolution(displacement) != -1)
            return false;
        assemble(displacement,assembleMatrix);
    }
    else // mixed formulation (displacement + pressure)
    {
        gsMultiPatch<T> displacement, pressure;
        constructSolution(solutionVector,displacement,pressure);
        if (checkSolution(displacement) != -1)
            return false;
        assemble(displacement,pressure,assembleMatrix);
    }
    return true;
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & displacement, bool assembleMatrix)
{

    if (assembleMatrix)
    {
        m_system.matrix().setZero();
        m_system.reserve(m_bases[0], m_options, 1);
    }
    m_system.rhs().setZero(Base::numDofs(),1);

    scaleDDoFs(m_options.getReal("DirichletAssembly"));

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNonLinearElasticity<T> visitor(*m_pde_ptr,displacement,assembleMatrix);
    Base::template push<gsVisitorNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    resetDDoFs();
    m_system.matrix().makeCompressed();
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & displacement, const gsMultiPatch<T> & pressure,
                                        bool assembleMatrix)
{
    if (assembleMatrix)
    {
        m_system.matrix().setZero();
        m_system.reserve(m_bases[0], m_options, 1);
    }
    m_system.rhs().setZero(Base::numDofs(),1);

    scaleDDoFs(m_options.getReal("DirichletAssembly"));

    // Compute volumetric integrals and write to the global linear system
    gsVisitorMixedNonLinearElasticity<T> visitor(*m_pde_ptr,displacement,pressure,assembleMatrix);
    Base::template push<gsVisitorMixedNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    resetDDoFs();
    m_system.matrix().makeCompressed();
}

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,result,unknowns);
}

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                                 gsMultiPatch<T>  & displacement, gsMultiPatch<T> & pressure) const
{
    GISMO_ENSURE(m_bases.size() == unsigned(m_dim) + 1, "Not a mixed formulation: can't construct pressure.");
    // construct displacement
    constructSolution(solVector,displacement);
    // construct pressure
    gsVector<index_t> unknowns(1);
    unknowns.at(0) = m_dim;
    Base::constructSolution(solVector,pressure,unknowns);
}

template <class T>
void gsElasticityAssembler<T>::constructPressure(const gsMatrix<T>& solVector, gsMultiPatch<T>& pressure) const
{
    gsVector<index_t> unknowns(1);
    unknowns.at(0) = m_dim;
    Base::constructSolution(solVector,pressure,unknowns);
}

template <class T>
void gsElasticityAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                       gsPiecewiseFunction<T> & result,
                                                       stress_type::type type) const
{
    // TODO: construct stresses for nonlinear and mixed elasticity
    result.clear();
    if (type == stress_type::all_2D)
        GISMO_ASSERT(m_dim == 2, "Invalid stress type for a 2D problem");
    if (type == stress_type::normal_3D || type == stress_type::shear_3D)
        GISMO_ASSERT(m_dim == 3, "Invalid stress type for a 3D problem");

    T E = m_options.getReal("YoungsModulus");
    T pr = m_options.getReal("PoissonsRatio");
    T lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
    T mu     = E / ( 2. * ( 1. + pr ) );

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(displacement,p,type,lambda,mu));
}

template <class T>
void gsElasticityAssembler<T>::deformGeometry(const gsMatrix<T> & solVector, gsMultiPatch<T> & domain) const
{
    GISMO_ASSERT(domain.domainDim() == m_dim,
                 "Wrong parametric dimension of a given domain: " + util::to_string(domain.domainDim()) +
                 ". Must be: " + util::to_string(m_dim) + ".\n");
    GISMO_ASSERT(domain.targetDim() == m_dim,
                 "Wrong dimension of a given domain: " + util::to_string(domain.targetDim()) +
                 ". Must be: " + util::to_string(m_dim) + ".\n");

    gsMultiPatch<T> solution;
    constructSolution(solVector,solution);

    GISMO_ASSERT(domain.nPatches() == solution.nPatches(),
                 "Wrong number of patches of a given domain: " + util::to_string(domain.nPatches()) +
                 ". Must be: " + util::to_string(solution.nPatches()) + ".\n");

    for (size_t p = 0; p < solution.nPatches(); ++p)
    {
        GISMO_ASSERT(domain.patch(p).coefsSize() == solution.patch(p).coefsSize(),
                     "Wrong number of control points in patch " + util::to_string(p) +
                     " of a given domain: " + util::to_string(domain.patch(p).coefsSize()) +
                     ". Must be: " + util::to_string(solution.patch(p).coefsSize()) + ".\n");

        domain.patch(p).coefs() += solution.patch(p).coefs();
    }
}

template <class T>
void gsElasticityAssembler<T>::setDirichletDofs(size_t patch, boxSide side, const gsMatrix<T> & ddofs)
{
    bool dirBcExists = false;
    typename gsBoundaryConditions<T>::const_iterator it = m_pde_ptr->bc().dirichletBegin();
    while (!dirBcExists && it != m_pde_ptr->bc().dirichletEnd())
    {
        if (it->patch() == patch && it->side() == side)
            dirBcExists = true;
        ++it;
    }
    GISMO_ASSERT(dirBcExists,"Side " + util::to_string(side) + " of patch " + util::to_string(patch)
                             + " does not belong to the Dirichlet boundary\n");

    gsMatrix<unsigned> localBIndices = m_bases[0][patch].boundary(side);
    GISMO_ASSERT(localBIndices.rows() == ddofs.rows() && m_dim == ddofs.cols(),
                 "Wrong size of a given matrix with Dirichlet DoFs: " + util::to_string(ddofs.rows()) +
                 " x " + util::to_string(ddofs.cols()) + ". Must be:" + util::to_string(localBIndices.rows()) +
                 " x " + util::to_string(m_dim) + ".\n");

    for (short_t d = 0; d < m_dim; ++d )
    {
        gsMatrix<unsigned> globalIndices;
        m_system.mapColIndices(localBIndices, patch, globalIndices, d);

        for (index_t i = 0; i < globalIndices.rows(); ++i)
            m_ddof[d](m_system.colMapper(d).global_to_bindex(globalIndices(i,0)),0) = ddofs(i,d);
    }
}

template <class T>
index_t gsElasticityAssembler<T>::checkSolution(const gsMultiPatch<T> & solution) const
{
    index_t corruptedPatch = -1;
    gsMapData<T> mdG, mdU;
    mdG.flags = NEED_DERIV;
    mdU.flags = NEED_DERIV;
    gsMatrix<T> points;

    for (size_t p = 0; p < solution.nPatches(); ++p)
    {
        gsQuadRule<T> quRule = gsQuadrature::get(m_bases[0][p],m_options);

        typename gsBasis<T>::domainIter domIt = m_bases[0][p].makeDomainIterator(boundary::none);
        for (; domIt->good(); domIt->next())
        {
            genSamplingPoints(domIt->lowerCorner(),domIt->upperCorner(),quRule,points);
            mdG.points = points;
            mdU.points = points;
            m_pde_ptr->domain().patch(p).computeMap(mdG);
            solution.patch(p).computeMap(mdU);
            for (index_t q = 0; q < points.cols(); ++q)
            {
                gsMatrix<T> physDispJac = mdU.jacobian(q)*(mdG.jacobian(q).cramerInverse());
                gsMatrix<T> F = gsMatrix<T>::Identity(m_dim,m_dim) + physDispJac;
                if (F.determinant() <= 0)
                    return p;
            }
        }
    }
    return corruptedPatch;
}

template <class T>
T gsElasticityAssembler<T>::solutionJacRatio(const gsMultiPatch<T> & solution) const
{
    std::vector<T> maxs, mins;
    gsMapData<T> mdG, mdU;
    mdG.flags = NEED_DERIV;
    mdU.flags = NEED_DERIV;
    gsMatrix<T> points;

    for (size_t p = 0; p < solution.nPatches(); ++p)
    {
        gsQuadRule<T> quRule = gsQuadrature::get(m_bases[0][p],m_options);

        typename gsBasis<T>::domainIter domIt = m_bases[0][p].makeDomainIterator(boundary::none);
        for (; domIt->good(); domIt->next())
        {
            genSamplingPoints(domIt->lowerCorner(),domIt->upperCorner(),quRule,points);
            mdG.points = points;
            mdU.points = points;
            m_pde_ptr->domain().patch(p).computeMap(mdG);
            solution.patch(p).computeMap(mdU);

            T minJ = (gsMatrix<T>::Identity(m_dim,m_dim) + mdU.jacobian(0)*(mdG.jacobian(0).cramerInverse())).determinant();
            T maxJ = minJ;
            for (int q = 1; q < points.cols(); ++q)
            {
                T J = (gsMatrix<T>::Identity(m_dim,m_dim) + mdU.jacobian(q)*(mdG.jacobian(q).cramerInverse())).determinant();
                if (J > maxJ)
                    maxJ = J;
                if (J < minJ)
                    minJ = J;
            }

            maxs.push_back(maxJ);
            mins.push_back(minJ);
        }
    }

    return *(std::min_element(mins.begin(),mins.end())) / *(std::max_element(maxs.begin(),maxs.end()));
}

template <class T>
void gsElasticityAssembler<T>::setDirichletAssemblyScaling(T factor)
{
    m_options.setReal("DirichletAssembly",factor);
}

template <class T>
void gsElasticityAssembler<T>::setDirichletConstructionScaling(T factor)
{
    m_options.setReal("DirichletConstruction",factor);
}

template <class T>
void gsElasticityAssembler<T>::setForceScaling(T factor)
{
    m_options.setReal("ForceScaling",factor);
}

template <class T>
void gsElasticityAssembler<T>::scaleDDoFs(T factor)
{
    if (saved_ddof.empty())
        saved_ddof = m_ddof;
    for (unsigned d = 0; d < m_ddof.size(); ++d)
        m_ddof[d] = saved_ddof[d]*factor;
}

template <class T>
void gsElasticityAssembler<T>::resetDDoFs()
{
    if (!saved_ddof.empty())
        m_ddof = saved_ddof;
}

template <class T>
void genSamplingPoints(const gsVector<T> & lower, const gsVector<T> & upper,
                       const gsQuadRule<T> & quRule, gsMatrix<T> & points)
{
    gsMatrix<T> quadPoints;
    gsVector<T> tempVector; // temporary argument for the gsQuadrule::mapTo function
    quRule.mapTo(lower,upper,quadPoints,tempVector);

    gsVector<unsigned> nPoints(quadPoints.rows());
    for (index_t d = 0; d < quadPoints.rows(); ++d)
        nPoints.at(d) = 2;
    gsMatrix<T> corners = gsPointGrid(lower,upper,nPoints);

    points.resize(quadPoints.rows(),quadPoints.cols()+corners.cols());
    points << quadPoints,corners;
}

}// namespace gismo ends
