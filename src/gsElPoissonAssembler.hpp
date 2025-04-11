/** @file gsElPoissonAssembler.hpp

    @brief Provides stiffness matrix for Poisson's equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/src/gsElPoissonAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsElasticity/src/gsVisitorElPoisson.h>

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
    opt.addReal("LocalStiff","Stiffening degree for the Jacobian-based local stiffening",0.);
    return opt;
}

template <class T>
void gsElPoissonAssembler<T>::refresh()
{
    std::vector<gsDofMapper> m_dofMappers(m_bases.size());
    m_bases[0].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                         iFace::glue,m_pde_ptr->bc(),m_dofMappers[0],0,true);

    m_system = gsSparseSystem<T>(m_dofMappers[0]);
    m_system.reserve(m_bases[0], m_options, m_pde_ptr->numRhs());
    Base::computeDirichletDofs(0);
}

template<class T>
void gsElPoissonAssembler<T>::assemble(bool saveEliminationMatrix)
{
    m_system.matrix().setZero();
    m_system.reserve(m_bases[0], m_options, m_pde_ptr->numRhs());
    m_system.rhs().setZero(Base::numDofs(),m_pde_ptr->numRhs());

    if (saveEliminationMatrix)
    {
        eliminationMatrix.resize(Base::numDofs(),Base::numFixedDofs());
        eliminationMatrix.setZero();
        eliminationMatrix.reservePerColumn(m_system.numColNz(m_bases[0],m_options));
    }

    gsVisitorElPoisson<T> visitor(*m_pde_ptr, saveEliminationMatrix ? &eliminationMatrix : nullptr);
    Base::template push<gsVisitorElPoisson<T> >(visitor);

    m_system.matrix().makeCompressed();

    if (saveEliminationMatrix)
    {
        Base::rhsWithZeroDDofs = m_system.rhs();
        eliminationMatrix.makeCompressed();
    }
}

template <class T>
void gsElPoissonAssembler<T>::constructSolution(const gsMatrix<T> & solVector,
                                                const std::vector<gsMatrix<T> > & fixedDoFs,
                                                gsMultiPatch<T> & result) const
{
    Base::constructSolution(solVector,fixedDoFs,result,gsVector<index_t>::Zero(1));
}

}// namespace gismo ends
