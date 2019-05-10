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
    // Always concieved as a meaningful class, now gsPde is just a container for
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

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

template <class T>
void gsElasticityAssembler<T>::refresh()
{
    GISMO_ENSURE(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ENSURE(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_dim);
    for (short_t d = 0; d < m_dim; d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_dim);
    m_system = gsSparseSystem<T>(m_dofMappers, dims);

    m_options.setReal("bdO",m_dim*(1+m_options.getReal("bdO"))-1);
    m_system.reserve(m_bases[0], m_options, 1);

    for (short_t d = 0; d < m_dim; ++d)
        Base::computeDirichletDofs(d);
}

template<class T>
void gsElasticityAssembler<T>::assemble(bool assembleMatrix)
{
    m_system.matrix().setZero();
    m_system.reserve(m_bases[0], m_options, 1);
    m_system.rhs().setZero(Base::numDofs(),1);

    if ( this->numDofs() == 0 )
    {
        gsWarn << "No internal DOFs. Computed Dirichlet boundary only.\n";
        return;
    }

    // Compute volumetric integrals and write to the global linear system
    gsVisitorLinearElasticity<T> visitor(*m_pde_ptr,assembleMatrix);
    Base::template push<gsVisitorLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed, bool assembleMatrix)
{
    if (assembleMatrix)
    {
        m_system.matrix().setZero();
        m_system.reserve(m_bases[0], m_options, 1);
    }
    m_system.rhs().setZero(Base::numDofs(),1);

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNonLinearElasticity<T> visitor(*m_pde_ptr,deformed,assembleMatrix);
    Base::template push<gsVisitorNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result, int unk) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,result,unknowns);
}

template <class T>
void gsElasticityAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                       gsPiecewiseFunction<T> & result,
                                                       stress_type::type type) const
{
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
    if (m_dim == 2)
        Base::constructSolution(solVector,solution,gsVector<index_t>::vec(0,1));
    else
        Base::constructSolution(solVector,solution,gsVector<index_t>::vec(0,1,2));

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
