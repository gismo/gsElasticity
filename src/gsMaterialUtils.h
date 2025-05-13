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
    /// @brief Deformation gradient
    F       = 00,
    /// @brief Deformation gradient in Voigt notation
    F_voigt = 01,
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

// Converts stresses into Voigt notation
// The stresses are assumed to be stored in a matrix of size d*d x n,
// where d is the dimension and n is the number evaluations.
template <class T>
void calculate_voigt_stress(const gsMatrix<T>& tmp, short_t d, gsMatrix<T>& result)
{
    GISMO_ASSERT(tmp.rows() == d*d, "Invalid size of stress matrix");
    result.resize(d*(d+1)/2, tmp.cols());
    for (index_t k = 0; k < tmp.cols(); k++)
    {
        gsMatrix<T> S = tmp.reshapeCol(k, d, d); // in tensor notation
        gsVector<T> S_voigt;
        voigtStress(S_voigt, S); // Convert to Voigt notation
        result.col(k) = S_voigt;
    }
}

template <class T>
void calculate_voigt_strain(const gsMatrix<T>& tmp, short_t d, gsMatrix<T>& result)
{
    GISMO_ASSERT(tmp.rows() == d*d, "Invalid size of stress matrix");
    result.resize(d*(d+1)/2, tmp.cols());
    for (index_t k = 0; k < tmp.cols(); k++)
    {
        gsMatrix<T> E = tmp.reshapeCol(k, d, d); // in tensor notation
        gsVector<T> E_voigt;
        voigtStrain(E_voigt, E); // Convert to Voigt notation
        result.col(k) = E_voigt;
    }
}

} // namespace
