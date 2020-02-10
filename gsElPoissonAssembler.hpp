/** @file gsElMassAssembler.hpp

    @brief Provides stiffness matrix for Poisson's equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElPoissonAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsElasticity/gsVisitorElPoisson.h>

namespace gismo
{

template<class T>
gsElPoissonAssembler<T>::gsElPoissonAssembler(const gsMultiPatch<T> & patches,
                                              const gsMultiBasis<T> & basis,
                                              const gsBoundaryConditions<T> & bconditions,
                                              const gsFunction<T> & body_force)
{
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    m_bases.push_back(basis);
    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsElPoissonAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("LocalStiff","Local stiffening",0.);
    return opt;
}

template <class T>
void gsElPoissonAssembler<T>::refresh()
{
    std::vector<gsDofMapper> m_dofMappers(m_bases.size());
    m_bases[0].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                         iFace::glue,m_pde_ptr->bc(),m_dofMappers[0],0,true);

    m_system = gsSparseSystem<T>(m_dofMappers[0]);
    m_system.reserve(m_bases[0], m_options, 1);
    Base::computeDirichletDofs(0);
}

template<class T>
void gsElPoissonAssembler<T>::assemble(bool assembleMatrix)
{
    m_system.matrix().setZero();
    m_system.reserve(m_bases[0], m_options, 1);
    m_system.rhs().setZero(Base::numDofs(),1);

    gsVisitorElPoisson<T> visitor(*m_pde_ptr);
    Base::template push<gsVisitorElPoisson<T> >(visitor);

    m_system.matrix().makeCompressed();
}

template <class T>
void gsElPoissonAssembler<T>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & result) const
{
    gsAssembler<T>::constructSolution(solVector,result);
}

}// namespace gismo ends
