/** @file gsNsAssembler.hpp

    @brief Implementation of gsNsAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsNsAssembler.h>

#include <gsPde/gsPoissonPde.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorStokes.h>
#include <gsElasticity/gsVisitorNavierStokes.h>

namespace gismo
{

template<class T>
gsNsAssembler<T>::gsNsAssembler(gsMultiPatch<T> const & patches,
                                gsMultiBasis<T> const & basisVel,
                                gsMultiBasis<T> const & basisPres,
                                gsBoundaryConditions<T> const & bconditions,
                                const gsFunction<T> & body_force)
{
    // comment: same problems as in gsElasticityAssembler
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    // same as above
    m_dim = body_force.targetDim();
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basisVel);
    m_bases.push_back(basisPres);

    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsNsAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("Viscosity","Kinematic viscosity of the fluid",0.001);
    opt.addReal("DirichletConstruction","Dirichlet BC scaling parameter for solution construction",1.);
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addReal("DirichletAssembly","Dirichlet BC scaling parameter for assembly",1.);
    opt.addSwitch("SUPG","Use SUPG stabilaztion",true);
    opt.addInt("Iteration","Type of the linear iteration used to solve the nonlinear problem",iteration_type::newton);
    return opt;
}

template <class T>
void gsNsAssembler<T>::refresh()
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
void gsNsAssembler<T>::assemble()
{
    m_system.matrix().setZero();
    m_system.reserve(m_bases[0], m_options, 1);
    m_system.rhs().setZero(Base::numDofs(),1);

    gsVisitorStokes<T> visitor(*m_pde_ptr);
    Base::template push<gsVisitorStokes<T> >(visitor);

    m_system.matrix().makeCompressed();
}

template <class T>
bool gsNsAssembler<T>::assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix)
{
    gsMultiPatch<T> velocity, pressure;
    constructSolution(solutionVector,velocity,pressure);

    if (assembleMatrix)
    {
        m_system.matrix().setZero();
        m_system.reserve(m_bases[0], m_options, 1);
    }
    m_system.rhs().setZero(Base::numDofs(),1);

    Base::scaleDDoFs(m_options.getReal("DirichletAssembly"));

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNavierStokes<T> visitor(*m_pde_ptr,velocity,pressure,assembleMatrix);
    Base::template push<gsVisitorNavierStokes<T> >(visitor);

    Base::resetDDoFs();
    m_system.matrix().makeCompressed();

    return true;
}

template <class T>
void gsNsAssembler<T>::constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& velocity) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,velocity,unknowns);
}

template <class T>
void gsNsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                         gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure) const
{
    // construct displacement
    constructSolution(solVector,velocity);
    // construct pressure
    constructPressure(solVector,pressure);
}

template <class T>
void gsNsAssembler<T>::constructPressure(const gsMatrix<T>& solVector, gsMultiPatch<T>& pressure) const
{
    gsVector<index_t> unknowns(1);
    unknowns.at(0) = m_dim;
    Base::constructSolution(solVector,pressure,unknowns);
}

} // namespace ends
