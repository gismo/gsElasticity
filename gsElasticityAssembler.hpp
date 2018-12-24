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

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
#include <gsElasticity/gsVisitorNonLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
                                                gsMultiBasis<T> const & basis,
                                                gsBoundaryConditions<T> const & bconditions,
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
    for (int d = 0; d < m_dim; ++d)
        m_bases.push_back(basis);

    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("Density","Density of the material",1.);
    opt.addReal("TimeFactor","Time factor for the time-dependent forces",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",0);
    return opt;
}

template <class T>
void gsElasticityAssembler<T>::refresh()
{
    GISMO_ASSERT(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_dim);
    for (index_t d = 0; d < m_dim; d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_dim);
    m_system = gsSparseSystem<T>(m_dofMappers, dims);

    m_system.reserve(m_bases[0], m_options, 1);

    for (int d = 0; d < m_dim; ++d)
        Base::computeDirichletDofs(d);
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{
    m_system.matrix().setZero();
    m_system.rhs().setZero();

    if ( this->numDofs() == 0 )
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n";
        return;
    }

    // Compute volumetric integrals and write to the global linear system
    Base::template push<gsVisitorLinearElasticity<T> >();
    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
    m_system.matrix().setZero();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNonLinearElasticity<T> visitor(*m_pde_ptr,deformed);
    Base::template push<gsVisitorNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
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

    for (index_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(displacement,p,type,lambda,mu));
}

template <class T>
void gsElasticityAssembler<T>::deformGeometry(const gsMatrix<T> & solVector, gsMultiPatch<T> & domain)
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

    for (index_t p = 0; p < solution.nPatches(); ++p)
    {
        GISMO_ASSERT(domain.patch(p).coefsSize() == solution.patch(p).coefsSize(),
                     "Wrong number of control points in patch " + util::to_string(p) +
                     " of a given domain: " + util::to_string(domain.patch(p).coefsSize()) +
                     ". Must be: " + util::to_string(solution.patch(p).coefsSize()) + ".\n");

        domain.patch(p).coefs() += solution.patch(p).coefs();
    }
}

template <class T>
void gsElasticityAssembler<T>::setDirichletDofs(index_t patch, boxSide side, const gsMatrix<T> & ddofs)
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

    for (index_t d = 0; d < m_dim; ++d )
    {
        gsMatrix<unsigned> globalIndices;
        m_system.mapColIndices(localBIndices, patch, globalIndices, d);

        for (index_t i = 0; i < globalIndices.rows(); ++i)
            m_ddof[d](m_system.colMapper(d).global_to_bindex(globalIndices(i,0)),0) = ddofs(i,d);
    }
}

}// namespace gismo ends
