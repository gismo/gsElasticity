/** @file gsElasticityAssembler.h

    @brief Provides linear and nonlinear elasticity systems for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

namespace gismo
{

/** @brief Assembles stiffness and mass matrices and right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. Matrices and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
*/
template <class T>
class gsElasticityAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

    /// @brief Constructor of the assembler object.
    gsElasticityAssembler(  gsMultiPatch<T> const & patches,
                            gsMultiBasis<T> const & bases,
                            gsBoundaryConditions<T> const & bconditions,
                            const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the stiffness matrix
    virtual void assemble();

    /// @brief Construct patchwise stress-function for visualization.
    void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                 gsPiecewiseFunction<T>& result,
                                 stress_type::type type = stress_type::von_mises) const;

protected:

    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
	index_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
