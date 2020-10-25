/** @file gsThermoAssembler.hpp

    @brief Provides assembler implementation for the gsThermoAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        D. Fusseder  (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsThermoAssembler.h>

#include <gsElasticity/gsVisitorThermo.h>
#include <gsElasticity/gsVisitorThermoBoundary.h>

namespace gismo
{

template <class T>
gsThermoAssembler<T>::gsThermoAssembler(const gsMultiPatch<T> & patches,
                                        const gsMultiBasis<T> & bases,
                                        const gsBoundaryConditions<T> & b_conditions,
                                        const gsFunction<T> & body_force,
                                        const gsFunctionSet<T> & temperature_field)
    :gsElasticityAssembler<T>(patches,bases,b_conditions,body_force),
     m_temperatureField(temperature_field),
     assembledElasticity(false)
{
    m_options.addReal("InitTemp","Initial temperature of the object",20.);
    m_options.addReal("ThExpCoef","Coefficient of thermal expansion of the material",20.);
    m_options.addSwitch("ParamTemp","Yes if the temperature field is parametric",true);

    findNonDirichletSides();
}

template <class T>
void gsThermoAssembler<T>::assemble(bool saveEliminationMatrix)
{
    gsElasticityAssembler<T>::assemble();
    elastRhs = gsAssembler<T>::m_system.rhs();
    assembledElasticity = true;

    assembleThermo();
}

template <class T>
void gsThermoAssembler<T>::findNonDirichletSides()
{
    for (std::vector< patchSide >::iterator side = m_pde_ptr->domain().bBegin(); side != m_pde_ptr->domain().bEnd(); ++side)
    {
        std::pair<index_t,boxSide> temp(side->patch,side->index());

        typename gsBoundaryConditions<T>::const_iterator it = m_pde_ptr->bc().dirichletBegin();
        for ( ; it != m_pde_ptr->bc().dirichletEnd(); ++it )
            if (temp.first == it->patch() && temp.second == it->side())
                goto exitLabel;

        nonDirichletSides.push_back(temp);
    exitLabel:;
    }
}

template <class T>
void gsThermoAssembler<T>::assembleThermo()
{
    GISMO_ENSURE(assembledElasticity, "gsElThermoAssembler::assemble() hasn't been called!");
    gsAssembler<T>::m_system.rhs().setZero();

    gsVisitorThermo<T> visitor(m_temperatureField);
    gsAssembler<T>::template push<gsVisitorThermo<T> >(visitor);

    for (std::vector<std::pair<index_t,boxSide> >::const_iterator it = nonDirichletSides.begin();
         it != nonDirichletSides.end(); ++it)
    {
        gsVisitorThermoBoundary<T> bVisitor(it->second,m_temperatureField);
        gsAssembler<T>::template apply<gsVisitorThermoBoundary<T> >(bVisitor,it->first,it->second);
    }
    gsAssembler<T>::m_system.rhs() += elastRhs;
}

} // namespace ends
