/** @file gsMixedElasticityAssembler.hpp

    @brief Provides linear and nonlinear mixed elasticity systems for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsMixedElasticityAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorMixedLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

namespace gismo
{

template<class T>
gsMixedElasticityAssembler<T>::gsMixedElasticityAssembler(gsMultiPatch<T> const & patches,
                                                          gsMultiBasis<T> const & basisDisplacement,
                                                          gsMultiBasis<T> const & basisPressure,
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
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basisDisplacement);
    m_bases.push_back(basisPressure);

    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsMixedElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::saint_venant_kirchhoff);
    return opt;
}

template <class T>
void gsMixedElasticityAssembler<T>::refresh()
{
    GISMO_ENSURE(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ENSURE(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_dim+1); // displacement + pressure
    for (short_t d = 0; d < m_dim; d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);
    m_bases.back().getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers.back(),m_dim,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_dim+1);
    m_system = gsSparseSystem<T>(m_dofMappers, dims);

    m_options.setReal("bdO",(m_dim+1)*(1+m_options.getReal("bdO"))-1);
    m_system.reserve(m_bases[0], m_options, 1);

    for (short_t d = 0; d < m_dim+1; ++d)
        Base::computeDirichletDofs(d);
}

template<class T>
void gsMixedElasticityAssembler<T>::assemble()
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
    Base::template push<gsVisitorMixedLinearElasticity<T> >();
    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template <class T>
void gsMixedElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>  & displacement, gsMultiPatch<T> & pressure) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,displacement,unknowns);
    Base::constructSolution(solVector,pressure,m_dim);
}

}// namespace gismo ends
