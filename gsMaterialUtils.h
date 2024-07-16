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
    /// @brief Positive part of the Strain
    Eplus   = 11,
    /// @brief Negative part of the Strain
    Eneg    = 12,
    /// @brief Stress
    S       = 2,
    /// @brief Material matrix
    C       = 3,
    /// @brief Energy
    Psi     = 4,
};

} // namespace
