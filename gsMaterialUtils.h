/** @file gsMaterialEval.h

    @brief Utilities for gsMaterial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

namespace gismo
{

enum class gsMaterialOutput : short_t
{
    /// @brief Strain
    E       = 10,
    /// @brief Strain in Voigt notation
    E_voigt = 11,
    /// @brief Stress
    S       = 20,
    /// @brief Stress in Voigt notation
    S_voigt = 21,
    /// @brief Material matrix in Voigt notation
    C       = 30,
    /// @brief Material matrix in Voigt notation undegraded
    C_pos   = 31,
    /// @brief Energy
    Psi     = 40,
    /// @brief Energy
    Psi_pos = 41,
};

} // namespace
