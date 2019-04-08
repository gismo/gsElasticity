/** @file gsElThermoAssembler.h

    @brief Provides a thermal expansion solver for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A.Shamanskiy
*/
#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

/** @brief Assembles stiffness and mass matrices and right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. Matrices and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
*/

// TODO: efficient RHS reassembly

template <class T>
class gsElThermoAssembler : public gsElasticityAssembler<T>
{
public:
    typedef gsElasticityAssembler<T> Base;

    /// @brief Constructor of the assembler object.
    gsElThermoAssembler(const gsMultiPatch<T> & patches,
                        const gsMultiBasis<T> & bases,
                        const gsBoundaryConditions<T> & b_conditions,
                        const gsFunction<T> & body_force,
                        const gsFunctionSet<T> & temperature_field);

    /// @brief Assembles the stiffness matrix and the RHS
    void assemble();

protected:
    /// @brief Marks all non-Dirichlet sides for assembly of the boundary thermal stresses
    void findNonDirichletSides();

    /// @brief Assembles the thermal expanstion contribution to the stiffness matrix and the RHS
    void assembleThermo(const gsFunctionSet<T> & temperatureField);

protected:
    const gsFunctionSet<T> & m_temperatureField;
    bool assembledElasticity;
    std::vector<std::pair<int,boxSide> > nonDirichletSides;

    using Base::m_pde_ptr;
    using Base::m_options;

}; // class definition ends
} // namespace ends
