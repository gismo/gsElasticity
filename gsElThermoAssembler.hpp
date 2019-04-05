/** @file gsElThermoAssembler.hpp

    @brief Provides assembler implementation for the gsElThermoAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsElasticity/gsElThermoAssembler.h>

#include <gsElasticity/gsVisitorElThermo.h>
#include <gsElasticity/gsVisitorElThermoBoundary.h>

namespace gismo
{

template <class T>
gsElThermoAssembler<T>::gsElThermoAssembler(const gsMultiPatch<T> & patches,
                                            const gsMultiBasis<T> & bases,
                                            const gsBoundaryConditions<T> & b_conditions,
                                            const gsFunction<T> & body_force,
                                            const gsFunctionSet<T> & heat_field)
    :gsElasticityAssembler<T>(patches,bases,b_conditions,body_force),
     m_heatField(heat_field),
     assembledElasticity(false)
{
    m_options.addReal("InitTemp","Initial temperature of the object",20.);
    m_options.addReal("ThExpCoef","Coefficient of thermal expansion of the material",20.);

    findNonDirichletSides();
}

template <class T>
void gsElThermoAssembler<T>::assemble()
{
    gsElasticityAssembler<T>::assemble();
    assembledElasticity = true;

    assembleThermo(m_heatField);
}

template <class T>
void gsElThermoAssembler<T>::findNonDirichletSides()
{
    for (std::vector< patchSide >::iterator side = m_pde_ptr->domain().bBegin(); side != m_pde_ptr->domain().bEnd(); ++side)
    {
        std::pair<int,int> temp(side->patch,side->index());

        typename gsBoundaryConditions<T>::const_iterator it = m_pde_ptr->bc().dirichletBegin();
        for ( ; it != m_pde_ptr->bc().dirichletEnd(); ++it )
            if (temp.first == it->patch() && temp.second == it->side())
                goto exitLabel;

        nonDirichletSides.push_back(temp);
        exitLabel:;
    }
}

template <class T>
void gsElThermoAssembler<T>::assembleThermo(const gsFunctionSet<T> & heatField)
{
    GISMO_ENSURE(assembledElasticity, "gsElThermoAssembler::assemble() hasn't been called!");
    /*m_rhsExtra = m_rhs;

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
    }*/
}

} // namespace ends
