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

#include <gsElasticity/gsBasePde.h>
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
    typename gsPde<T>::Ptr pde( new gsBasePde<T>(patches,bconditions,rightHandSides) );

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
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addInt("Assembly","Type of the linear system to assemble",ns_assembly::newton_update);
    return opt;
}

template <class T>
void gsNsAssembler<T>::reserve()
{
    // Pick up values from options
    const T bdA       = m_options.getReal("bdA");
    const index_t bdB = m_options.getInt("bdB");
    const T bdO       = m_options.getReal("bdO");

    index_t deg = 0;
    for (index_t d = 0; d < m_bases[0][0].dim(); ++d )
        if (m_bases[0][0].degree(d) > deg)
            deg = m_bases[0][0].degree(d);

    // m_dim velocity*velocity blocks + 1 pressure*velocity block (additioanal factor 2 for subgrid element)
    index_t numElPerColumn = pow((bdA*deg+bdB),m_dim)*m_dim + pow((2*bdA*deg+bdB),m_dim);

    m_system.reserve(numElPerColumn*(1+bdO),1);
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

    m_system = gsSparseSystem<T>(m_dofMappers, gsVector<index_t>::Ones(m_bases.size()));
    reserve();

    for (unsigned d = 0; d < m_bases.size(); ++d)
        Base::computeDirichletDofs(d);
}

//--------------------- SYSTEM ASSEMBLY ----------------------------------//

template<class T>
void gsNsAssembler<T>::assemble(bool saveEliminationMatrix)
{
    GISMO_UNUSED(saveEliminationMatrix);
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    gsVisitorStokes<T> visitor(*m_pde_ptr);
    Base::template push<gsVisitorStokes<T> >(visitor);

    m_system.matrix().makeCompressed();
}

template <class T>
bool gsNsAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
                                const std::vector<gsMatrix<T> > & fixedDoFs)
{
    gsMultiPatch<T> velocity, pressure;
    constructSolution(solutionVector,fixedDoFs,velocity,pressure);
    assemble(velocity,pressure);

    return true;
}

template <class T>
void gsNsAssembler<T>::assemble(const gsMultiPatch<T> & velocity,
                                const gsMultiPatch<T> & pressure)
{
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    gsVisitorNavierStokes<T> visitor(*m_pde_ptr,velocity,pressure);
    Base::template push<gsVisitorNavierStokes<T> >(visitor);

    m_system.matrix().makeCompressed();
}

//--------------------- SOLUTION CONSTRUCTION ----------------------------------//

template <class T>
void gsNsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                         const std::vector<gsMatrix<T> > & fixedDoFs,
                                         gsMultiPatch<T>& velocity) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,fixedDoFs,velocity,unknowns);
}

template <class T>
void gsNsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                         const std::vector<gsMatrix<T> > & fixedDoFs,
                                         gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure) const
{
    // construct displacement
    constructSolution(solVector,fixedDoFs,velocity);
    // construct pressure
    constructPressure(solVector,fixedDoFs,pressure);
}

template <class T>
void gsNsAssembler<T>::constructPressure(const gsMatrix<T>& solVector,
                                         const std::vector<gsMatrix<T> > & fixedDoFs,
                                         gsMultiPatch<T>& pressure) const
{
    GISMO_UNUSED(fixedDoFs);
    gsVector<index_t> unknowns(1);
    unknowns.at(0) = m_dim;
    Base::constructSolution(solVector,m_ddof,pressure,unknowns);
}

//--------------------- SPECIALS ----------------------------------//

template <class T>
gsMatrix<T> gsNsAssembler<T>::computeForce(const gsMultiPatch<T> & velocity, const gsMultiPatch<T> & pressure,
                                           const std::vector<std::pair<index_t,boxSide> > & bdrySides, bool split) const
{
    // all temporary data structures
    gsMatrix<T> quNodes, pressureValues, physGradJac, sigma;
    gsVector<T> quWeights, normal;
    // NEED_MEASURE for integration
    // NEED_GRAD_TRANSFORM for velocity gradients transformation from parametric to physical domain
    gsMapData<T> mdGeo(NEED_MEASURE | NEED_GRAD_TRANSFORM);
    // NEED_DERIV for velocity gradients
    gsMapData<T> mdVel(NEED_DERIV);

    gsMatrix<T> force;
    force.setZero(m_dim, 2);
    const T viscosity = m_options.getReal("Viscosity");
    const T density = m_options.getReal("Density");

    // loop over bdry sides
    for (std::vector<std::pair<index_t,boxSide> >::const_iterator it = bdrySides.begin();
         it != bdrySides.end(); ++it)
    {
        // basis of the patch
        const gsBasis<T> & basis = m_bases[0][it->first];
        // setting quadrature rule for the boundary side
        gsGaussRule<T> bdQuRule(basis,1.0,1,it->second.direction());
        // loop over elements of the side
        typename gsBasis<T>::domainIter elem = basis.makeDomainIterator(it->second);
        for (; elem->good(); elem->next())
        {
            // mapping quadrature rule to the element
            bdQuRule.mapTo(elem->lowerCorner(),elem->upperCorner(),quNodes,quWeights);
            // evaluate geoemtry mapping at the quad points
            mdGeo.points = quNodes;
            m_pde_ptr->patches().patch(it->first).computeMap(mdGeo);
            // evaluate velocity at the quad points
            mdVel.points = quNodes;
            velocity.patch(it->first).computeMap(mdVel);
            // evaluate pressure at the quad points
            pressure.patch(it->first).eval_into(quNodes,pressureValues);

            // loop over quad points
            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                // transform gradients from parametric to physical
                physGradJac = mdVel.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
                // normal length is the local measure
                outerNormal(mdGeo,q,it->second,normal);
                // pressure contribution to the stress tensor
                sigma = pressureValues.at(q)*gsMatrix<T>::Identity(m_dim,m_dim);
                force.col(0) += quWeights[q] * sigma * normal;
                sigma = -1*density*viscosity*(physGradJac + physGradJac.transpose());
                force.col(1) += quWeights[q] * sigma * normal;
            }
        }
    }
    if (split)
        return force;
    else
        return force.col(0) + force.col(1);
}

} // namespace ends
