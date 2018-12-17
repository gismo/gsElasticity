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

#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsField.h>
#include <gsPde/gsPoissonPde.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
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
    opt.addReal("Density","Density of the material",8e3);
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
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{
    if ( this->numDofs() == 0 )
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n";
        return;
    }

    m_system.reserve(m_bases[0], m_options, 1);

    for (int d = 0; d < m_dim; ++d)
        Base::computeDirichletDofs(d);

    // Compute volumetric integrals and write to the global linear system
    Base::template push<gsVisitorLinearElasticity<T> >();
    // Compute surface integrals and write to the global rhs vector
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

}// namespace gismo ends
