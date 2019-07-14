/** @file gsElUtils.h

    @brief Provides several simple utility classes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>

namespace gismo
{

struct newton_verbosity
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2  /// output everything
    };
};

struct newton_save
{
    enum save
    {
        onlyFinal = 0,  /// save only the final solution
        firstAndLastPerIncStep = 1,  /// save only the first and the last displacement fields at each incremental loading step
        all = 2  /// save every intermediate displacement field
    };
};

/** @brief Specifies the type of stresses to compute
 *
 *         Currently, gsWriteParaview can only plot vector-valued functions with an output dimension up to three.
 *         Therefore it not possible to plot all stress components as components of a single vector-valued function.
*/
struct stress_type
{
    enum type
    {
        von_mises = 0,  /// compute only von Mises stress
        all_2D    = 1,  /// compute normal and shear stresses in 2D case (s11 s22 s12)
        normal_3D = 2,  /// compute normal stresses in 3D case (s11 s22 s33)
        shear_3D  = 3   /// compute shear stresses in 3D case (s12 s13 s23)
    };
};

/// @brief Specifies the material law to use
struct material_law
{
    enum type
    {
        saint_venant_kirchhoff = 0,  /// S = 2*mu*E + lambda*tr(E)*I
        neo_hooke_ln           = 1,  /// S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
        neo_hooke_2            = 2   /// S = lambda/2*(J^2-1)*C^-1 + mu*(I-C^-1)
    };
};

enum class elasticity_formulation { displacement, mixed_pressure };


enum class newton_status { converged, interrupted, working, bad_solution };

class gsProgressBar
{
public:
    gsProgressBar(index_t width = 25) : m_width(width) {}

    void display(double progress)
    {
        GISMO_ENSURE(progress >= 0. && progress <= 1.,"Invalid progress value! Must be between 0 and 1.");

        gsInfo << "[";
        for(index_t i = 0; i < m_width; i++)
            if(i < index_t(progress*m_width))
                gsInfo << "=";
            else if(i == index_t(progress*m_width))
                gsInfo << ">";
            else
                gsInfo << " ";
        gsInfo << "] " << progress*100 << " %\r";
        gsInfo.flush();

        if (abs(progress - 1.) < 1e-12)
            gsInfo << std::endl;
    }
protected:
    index_t m_width;
};

}
