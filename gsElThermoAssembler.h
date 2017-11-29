/** @file gsElThermoAssembler.h

   @brief Provides assembler for the thermo-elasticity equation.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

template <class T>
class gsElThermoAssembler : public gsElasticityAssembler<T>
{
public:

    gsElThermoAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & bases,
                        T youngModulus,
                        T poissonsRatio,
                        T thermalExpCoef,
                        T initTemp,
                        const gsBoundaryConditions<T> & bcInfo,
                        const gsFunction<T> & force,
                        const gsMultiPatch<T> & heatSolution,
                        dirichlet::strategy enforceStrategy = dirichlet::elimination,
                        dirichlet::values computeStrategy = dirichlet::l2Projection);

    /// Main assembling routine. Assemble elasticity contribution,
    /// as well as thermo contribution. The latter can be recomputed
    /// independently by "assembleThermo()",
    /// useful for time-dependent themo-elasticity
    void assemble();
    void setHeatSolution(const gsMultiPatch<T> & heatSolution);

protected:

    void findNonDirichletSides();
    void assembleThermo(const gsMultiPatch<T> & heatField);

protected:
    T m_thExpCoef;
    T m_initTemp;
    const gsMultiPatch<T> & m_heatSolution;
    std::vector<std::pair<int,int> > nonDirichletSides;

    using gsElasticityAssembler<T>::m_patches;
    using gsElasticityAssembler<T>::m_bConditions;
    using gsElasticityAssembler<T>::m_rhsExtra;
    using gsElasticityAssembler<T>::m_rhs;
    using gsElasticityAssembler<T>::m_lambda;
    using gsElasticityAssembler<T>::m_mu;

}; // class definition ends
} // namespace ends
