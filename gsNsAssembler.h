/** @file gsNsAssembler.h

    @brief Provides matrix and rhs assebmly for stationary and transient incompressible
    Stokes and Navier-Stokes equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>
#include <gsElasticity/gsBaseUtils.h>

namespace gismo
{

/** @brief TODO: write.
*/
template <class T>
class gsNsAssembler : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;

    /// @brief Constructor
    gsNsAssembler(const gsMultiPatch<T> & patches,
                  const gsMultiBasis<T> & basisVel,
                  const gsMultiBasis<T> & basisPres,
                  const gsBoundaryConditions<T> & bconditions,
                  const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembly of the linear system for the Stokes problem
    virtual void assemble();

    /// @ brief TODO: write
    virtual bool assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix = true) {return true;}

    /// @brief Construct velocity from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& velocity) const;

    /// @brief Construct velocity and pressure from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure) const;

    /// @ brief Construct pressure from computed solution vector
    virtual void constructPressure(const gsMatrix<T> & solVector, gsMultiPatch<T> & pressure) const;

    /// sets scaling of Dirichlet BC used for linear system assembly
    virtual void setDirichletAssemblyScaling(T factor) {}
    /// sets scaling of Dirichlet BC used for construction of the solution as a gsMultiPatch object
    virtual void setDirichletConstructionScaling(T factor) {}
    /// set scaling of the force loading (volume and surface loading)
    virtual void setForceScaling(T factor) {}

protected:

    /// Dimension of the problem
    /// parametric dim = physical dim = velocity dim
    short_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNsAssembler.hpp)
#endif
