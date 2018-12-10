/** @file gsElThermoAssembler.hpp

    @brief Provides assembler implementation for the gsElThermoAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsThermoElasticity/gsElThermoAssembler.h>

#include <gsThermoElasticity/gsVisitorElThermo.h>
#include <gsThermoElasticity/gsVisitorElThermoBoundary.h>

namespace gismo
{

template <class T>
gsElThermoAssembler<T>::gsElThermoAssembler(const gsMultiPatch<T> & patches,
                                            const gsMultiBasis<T> & bases,
                                            T youngModulus,
                                            T poissonsRatio,
                                            T thermalExpCoef,
                                            T initTemp,
                                            const gsBoundaryConditions<T> & bcInfo,
                                            const gsFunction<T> & force,
                                            //const gsMultiPatch<T> & heatSolution,
                                            const gsFunctionSet<T> & heatSolution,
                                            dirichlet::strategy enforceStrategy,
                                            dirichlet::values computeStrategy)
    :gsElasticityAssembler<T>(patches,bases,youngModulus,poissonsRatio,1.,bcInfo,force,computeStrategy,enforceStrategy),
     m_thExpCoef(thermalExpCoef), m_initTemp(initTemp),
     m_heatSolution(heatSolution)
{
    findNonDirichletSides();
}

template <class T>
void gsElThermoAssembler<T>::assemble()
{
    gsElasticityAssembler<T>::assemble();

    assembleThermo(m_heatSolution);
}

template <class T>
void gsElThermoAssembler<T>::findNonDirichletSides()
{
    for (std::vector< patchSide >::iterator side = m_patches.bBegin(); side != m_patches.bEnd(); ++side)
    {
        std::pair<int,int> temp(side->patch,side->index());

        for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.dirichletBegin();
              it != m_bConditions.dirichletEnd(); ++it )
        {
            if (temp.first == it->patch() && temp.second == it->side())
            {
                goto exitLabel;
            }
        }
        nonDirichletSides.push_back(temp);
        exitLabel:;
    }
}

template <class T>
//void gsElThermoAssembler<T>::assembleThermo(const gsMultiPatch<T> & heatField)
void gsElThermoAssembler<T>::assembleThermo(const gsFunctionSet<T> & heatField)
{
    GISMO_ASSERT(m_rhs.rows() != 0,
                 "Assemble() hasn't been called!");
    m_rhsExtra = m_rhs;

    for (index_t p = 0; p < m_patches.nPatches(); ++p)
    {
        gsVisitorElThermo<T> visitor(heatField.function(p),m_rhsExtra,
                                     m_lambda,m_mu,m_thExpCoef);
        this->apply(visitor,p);
    }

    for(std::vector<std::pair<int, int> >::iterator it = nonDirichletSides.begin(); it != nonDirichletSides.end(); ++it)
    {
        gsVisitorElThermoBoundary<T> bVisitor(heatField.function(it->first),it->second,m_rhsExtra,m_initTemp,
                                              m_lambda,m_mu,m_thExpCoef);
        this->apply(bVisitor,it->first,it->second);
    }
}

template <class T>
//void gsElThermoAssembler<T>::setHeatSolution(const gsMultiPatch<T> & heatSolution)
void gsElThermoAssembler<T>::setHeatSolution(const gsFunctionSet<T> & heatSolution)

{
    assembleThermo(heatSolution);
}



} // namespace ends
