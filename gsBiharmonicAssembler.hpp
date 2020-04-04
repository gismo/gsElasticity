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

#include <gsElasticity/gsBiharmonicAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsElasticity/gsVisitorBiharmonic.h>

namespace gismo
{
template<class T>
gsBiharmonicAssembler<T>::gsBiharmonicAssembler(const gsMultiPatch<T> & patches,
                                                const gsMultiBasis<T> & basis,
                                                const gsBoundaryConditions<T> & bconditions,
                                                const gsFunction<T> & body_force)
{
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    m_bases.push_back(basis);
    m_bases.push_back(basis);
    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsBiharmonicAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("LocalStiff","Stiffening degree for the Jacobian-based local stiffening",0.);
    return opt;
}

template <class T>
void gsBiharmonicAssembler<T>::refresh()
{
    std::vector<gsDofMapper> m_dofMappers(2);
    for (unsigned d = 0; d < 2; d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    m_system = gsSparseSystem<T>(m_dofMappers, gsVector<T>::Ones(2));
    reserve();
    Base::computeDirichletDofs(0);
    Base::computeDirichletDofs(1);
}

template <class T>
void gsBiharmonicAssembler<T>::reserve()
{
    // Pick up values from options
    const T bdA       = m_options.getReal("bdA");
    const index_t bdB = m_options.getInt("bdB");
    const T bdO       = m_options.getReal("bdO");

    index_t dim = m_bases[0][0].dim();
    index_t deg = 0;
    for (index_t d = 0; d < dim; ++d )
        if (m_bases[0][0].degree(d) > deg)
            deg = m_bases[0][0].degree(d);

    index_t numElPerColumn = pow((bdA*deg+bdB),dim)*2;
    m_system.reserve(numElPerColumn*(1+bdO),1);
}

template<class T>
void gsBiharmonicAssembler<T>::assemble()
{
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero(Base::numDofs(),m_pde_ptr->numRhs());

    gsVisitorBiharmonic<T> visitor(*m_pde_ptr);
    Base::template push<gsVisitorBiharmonic<T> >(visitor);

    m_system.matrix().makeCompressed();
}

template <class T>
void gsBiharmonicAssembler<T>::constructSolution(const gsMatrix<T> & solVector,
                                                const std::vector<gsMatrix<T> > & fixedDoFs,
                                                gsMultiPatch<T> & solutionMain, gsMultiPatch<T> & solutionAux) const
{
    Base::constructSolution(solVector,fixedDoFs,solutionMain,gsVector<index_t>::Zero(1));
    Base::constructSolution(solVector,fixedDoFs,solutionAux,gsVector<index_t>::Ones(1));
}

}// namespace gismo ends
