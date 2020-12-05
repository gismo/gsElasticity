/** @file gsMassAssembler.h

    @brief Provides mass matrix for elasticity systems in 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>

namespace gismo
{

/** @brief Assembles the mass matrix and right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. The matrix has a block structure associated with
           components of the displacement vector, each block corresponding to one component.
           Supports mixed displacement-pressure formulation.
*/

template <class T>
class gsMassAssembler : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;

    gsMassAssembler(const gsMultiPatch<T> & patches,
                      const gsMultiBasis<T> & basis,
                      const gsBoundaryConditions<T> & bconditions,
                      const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the mass matrix
    virtual void assemble(bool saveEliminationMatrix = false);

    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDDoFs)
    {assemble(); return true;}

protected:
    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
    short_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_options;
    using Base::m_system;
    using Base::m_ddof;
    using Base::eliminationMatrix;
};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElMassAssembler.hpp)
#endif
