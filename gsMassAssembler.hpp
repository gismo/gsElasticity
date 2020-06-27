/** @file gsElMassAssembler.hpp

    @brief Provides mass matrix for elasticity systems in 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsBasePde.h>

#include <gsPde/gsPoissonPde.h>
#include <gsElasticity/gsVisitorMass.h>

namespace gismo
{

template<class T>
gsMassAssembler<T>::gsMassAssembler(const gsMultiPatch<T> & patches,
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
    typename gsPde<T>::Ptr pde( new gsBasePde<T>(patches,bconditions,rightHandSides) );
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
gsOptionList gsMassAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("Density","Density of the material",1.);
    return opt;
}

template <class T>
void gsMassAssembler<T>::refresh()
{
    GISMO_ENSURE(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ENSURE(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_bases.size());
    for (unsigned d = 0; d < m_bases.size(); d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    m_system = gsSparseSystem<T>(m_dofMappers,gsVector<index_t>::Ones(m_bases.size()));

    m_options.setReal("bdO",m_bases.size()*(1+m_options.getReal("bdO"))-1);
    m_system.reserve(m_bases[0], m_options, 1);

    for (unsigned d = 0; d < m_bases.size(); ++d)
        Base::computeDirichletDofs(d);
}

template<class T>
void gsMassAssembler<T>::assemble(bool saveEliminationMatrix)
{
    // allocate space for the linear system
    m_system.matrix().setZero();
    m_system.reserve(m_bases[0], m_options, 1);
    m_system.rhs().setZero(Base::numDofs(),1);

    if (saveEliminationMatrix)
    {
        eliminationMatrix.resize(Base::numDofs(),Base::numFixedDofs());
        eliminationMatrix.setZero();
        eliminationMatrix.reservePerColumn(m_system.numColNz(m_bases[0],m_options));
    }

    gsVisitorMass<T> visitor(saveEliminationMatrix ? &eliminationMatrix : nullptr);
    Base::template push<gsVisitorMass<T> >(visitor);

    m_system.matrix().makeCompressed();

    if (saveEliminationMatrix)
    {
        Base::rhsWithZeroDDofs = m_system.rhs();
        eliminationMatrix.makeCompressed();
    }
}


}// namespace gismo ends
