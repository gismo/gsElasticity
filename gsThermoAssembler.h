/** @file gsThermoAssembler.h

    @brief Provides a thermal expansion solver for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        D. Fusseder  (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/
#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

/** @brief Assembles stiffness and mass matrices and right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. Matrices and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
*/

template <class T>
class gsThermoAssembler : public gsElasticityAssembler<T>
{
public:
    typedef gsElasticityAssembler<T> Base;

    /// @brief Constructor of the assembler object.
    gsThermoAssembler(const gsMultiPatch<T> & patches,
                      const gsMultiBasis<T> & bases,
                      const gsBoundaryConditions<T> & b_conditions,
                      const gsFunction<T> & body_force,
                      const gsFunctionSet<T> & temperature_field);

    /// @brief Assembles the stiffness matrix and the RHS
    virtual void assemble(bool saveEliminationMatrix = false);

    /// @brief Assembles the thermal expanstion contribution to the RHS
    void assembleThermo();

protected:
    /// @brief Marks all non-Dirichlet sides for assembly of the boundary thermal stresses
    void findNonDirichletSides();

protected:
    const gsFunctionSet<T> & m_temperatureField;
    bool assembledElasticity;
    std::vector<std::pair<int,boxSide> > nonDirichletSides;
    //
    gsMatrix<T> elastRhs;

    using Base::m_pde_ptr;
    using Base::m_options;

}; // class definition ends
} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsThermoAssembler.hpp)
#endif
