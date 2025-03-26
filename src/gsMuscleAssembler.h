/** @file gsMuscleAssembler.h

    @brief Provides elasticity systems for muscle simulations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

/** @brief Assembler for incompressible nonlinear elasticity problem with a muscle material model.
 * The material model is based on the paper by M.H.Gfrerer and B.Simeon
    "Fiber-based modeling and simulation of skeletal muscles"
*/
template <class T>
class gsMuscleAssembler : public gsElasticityAssembler<T>
{
public:
    typedef gsElasticityAssembler<T> Base;

    /// @brief Constructor of mixed formulation (displacement + pressure)
    gsMuscleAssembler(const gsMultiPatch<T> & patches,
                      const gsMultiBasis<T> & basisDisp,
                      const gsMultiBasis<T> & basisPres,
                      const gsBoundaryConditions<T> & bconditions,
                      const gsFunction<T> & body_force,
                      const gsPiecewiseFunction<T> & tendonMuscleDistribution,
                      const gsVector<T> & fiberDirection);

    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// Assembles the tangential linear system for Newton's method given the current solution
    /// in the form of free and fixed/Dirichelt degrees of freedom.
    /// Checks if the current solution is valid (Newton's solver can exit safely if invalid).
    using Base::assemble;
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs);

    //--------------------- SPECIALS ----------------------------------//

    /// @brief Construct Cauchy stresses for evaluation or visualization
    using Base::constructCauchyStresses;
    virtual void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                         const gsMultiPatch<T> & pressure,
                                         gsPiecewiseFunction<T> & result,
                                         stress_components::components component = stress_components::von_mises) const;

protected:
    using Base::m_options;
    using Base::m_pde_ptr;
    using Base::m_system;
    gsPiecewiseFunction<T> const & muscleTendon;
    gsVector<T> const & fiberDir;
};


} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMuscleAssembler.hpp)
#endif
