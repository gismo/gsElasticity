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

enum gsNeedMaterial
{
    NEED_MATERIAL_F        = 1U << 0, ///< Value of the object
    NEED_MATERIAL_E        = 1U << 1, ///< Gradient of the object
};

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

/**
 * @brief      Material data container
 *             This class contains deformation gradients and strains
 * @tparam     T     Real type
 * @ingroup    Elasticity
 */
template<class T>
class gsMaterialData
{

public:

    void reset()
    {
        deformationGradient.resize(0,0);
        strain.resize(0,0);
        parameters.clear();
        density.resize(0,0);
        dim = 0;
        size = 0;
        patch = 0;
        flags = 0;
    }

    void resizeParameters(index_t n)
    {
        parameters.resize(n);
    }

    // typename gsMaterialBase<T>::function_ptr m_undeformed;
    // typename gsMaterialBase<T>::function_ptr m_deformed;

    mutable std::vector<gsMatrix<T>> parameters; ///< Material parameters for each patch
    mutable             gsMatrix<T>  density; ///< Density for each patch

    // mutable gsMatrix<T> parameters;
    // mutable gsMatrix<T> m_rhoMat;
    // mutable gsMatrix<T> m_jac_ori, m_jac_def;
    mutable gsMatrix<T> deformationGradient, strain;

    mutable short_t dim;
    mutable index_t size;
    mutable index_t patch;
    mutable unsigned flags;
};

} // namespace
