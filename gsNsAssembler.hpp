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
    opt.addReal("Density","Density of the fluid",1.);
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
    if (!ALE)
    {
        gsVisitorNavierStokes<T> visitor(*m_pde_ptr,velocity,pressure,assembleMatrix);
        Base::template push<gsVisitorNavierStokes<T> >(visitor);
    }
    else
    {
        gsVisitorNavierStokes<T> visitor(*m_pde_ptr,velocity,pressure,*aleDisp,assembleMatrix);
        Base::template push<gsVisitorNavierStokes<T> >(visitor);
    }

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

template <class T>
gsMatrix<T> gsNsAssembler<T>::computeForce(const gsMultiPatch<T> & velocity, const gsMultiPatch<T> & pressure,
                                           const std::vector<std::pair<index_t,boxSide> > & bdrySides) const
{
    gsMatrix<T> force;
    force.setZero(m_dim,1);
    const T viscosity = m_options.getReal("Viscosity");
    const T density = m_options.getReal("Density");

    // loop over bdry sides
    for (auto &it : bdrySides)
    {
        // basis of the patch
        const gsBasis<T> & basis = m_bases[0][it.first];
        // setting quadrature rule for the boundary side
        gsGaussRule<T> bdQuRule(basis,1.0,1,it.second.direction());
        // loop over elements of the side
        typename gsBasis<T>::domainIter elem = basis.makeDomainIterator(it.second);
        for (; elem->good(); elem->next())
        {
            // mapping quadrature rule to the element
            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            bdQuRule.mapTo(elem->lowerCorner(),elem->upperCorner(),quNodes,quWeights);
            // evaluate geoemtry mapping at the quad points
            // NEED_MEASURE for integration
            // NEED_GRAD_TRANSFORM for velocity gradients transformation from parametric to physical domain
            gsMapData<T> mdGeo(NEED_MEASURE | NEED_GRAD_TRANSFORM);
            mdGeo.points = quNodes;
            m_pde_ptr->patches().patch(it.first).computeMap(mdGeo);
            // evaluate velocity at the quad points
            // NEED_DERIV for velocity gradients
            gsMapData<T> mdVel(NEED_DERIV);
            mdVel.points = quNodes;
            velocity.patch(it.first).computeMap(mdVel);
            // evaluate pressure at the quad points
            gsMatrix<T> pressureValues;
            pressure.patch(it.first).eval_into(quNodes,pressureValues);

            // loop over quad points
            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                // transform gradients from parametric to physical
                gsMatrix<T> physGradJac = mdVel.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
                // normal length is the local measure
                gsVector<T> normal;
                outerNormal(mdGeo,q,it.second,normal);
                // stress tensor
                gsMatrix<T> sigma = pressureValues.at(q)*gsMatrix<T>::Identity(m_dim,m_dim) -
                                    density*viscosity*(physGradJac + physGradJac.transpose());
                force += quWeights[q] * sigma * normal;
            }
        }
    }
    return force;
}

template <class T>
gsMatrix<T> gsNsAssembler<T>::computeForceALE(const gsMultiPatch<T> & velocity, const gsMultiPatch<T> & pressure,
                                              const gsMultiPatch<T> & displacementALE,
                                              const std::vector<std::pair<index_t,boxSide> > & bdrySides) const
{
    gsMatrix<T> force;
    force.setZero(m_dim,1);
    const T viscosity = m_options.getReal("Viscosity");
    const T density = m_options.getReal("Density");

    // loop over bdry sides
    for (auto &it : bdrySides)
    {
        // basis of the patch
        const gsBasis<T> & basis = m_bases[0][it.first];
        // setting quadrature rule for the boundary side
        gsGaussRule<T> bdQuRule(basis,1.0,1,it.second.direction());
        // loop over elements of the side
        typename gsBasis<T>::domainIter elem = basis.makeDomainIterator(it.second);
        for (; elem->good(); elem->next())
        {
            // mapping quadrature rule to the element
            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            bdQuRule.mapTo(elem->lowerCorner(),elem->upperCorner(),quNodes,quWeights);
            // evaluate geoemtry mapping at the quad points
            // NEED_MEASURE for integration
            // NEED_GRAD_TRANSFORM for velocity gradients transformation from parametric to physical domain
            gsMapData<T> mdGeo(NEED_MEASURE | NEED_GRAD_TRANSFORM);
            mdGeo.points = quNodes;
            m_pde_ptr->patches().patch(it.first).computeMap(mdGeo);
            // evaluate velocity at the quad points
            // NEED_DERIV for velocity gradients
            gsMapData<T> mdVel(NEED_DERIV);
            mdVel.points = quNodes;
            velocity.patch(it.first).computeMap(mdVel);
            // evaluate pressure at the quad points
            gsMatrix<T> pressureValues;
            pressure.patch(it.first).eval_into(quNodes,pressureValues);
            // evaluate ALE mapping at the param points
            // NEED_DERIV for gradients
            gsMapData<T> mdALE(NEED_DERIV);
            mdALE.points = quNodes;
            displacementALE.patch(it.first).computeMap(mdALE);

            // loop over quad points
            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                // transform gradients from parametric to physical
                gsMatrix<T> physGradJac = mdVel.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
                // ALE jacobian (identity + physical displacement gradient)
                gsMatrix<T> physJacALE = gsMatrix<T>::Identity(m_dim,m_dim) +
                        mdALE.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
                // inverse ALE jacobian
                gsMatrix<T> invJacALE = physJacALE.cramerInverse();
                // normal length is the local measure
                gsVector<T> normal;
                outerNormal(mdGeo,q,it.second,normal);
                // stress tensor
                gsMatrix<T> sigma = pressureValues.at(q)*gsMatrix<T>::Identity(m_dim,m_dim) -
                                    density*viscosity*(physGradJac*invJacALE +
                                                       invJacALE.transpose()*physGradJac.transpose());
                // stress tensor pull back
                gsMatrix<T> sigmaALE = physJacALE.determinant()*sigma*(invJacALE.transpose());
                force += quWeights[q] * sigmaALE * normal;
            }
        }
    }
    return force;
}

} // namespace ends
