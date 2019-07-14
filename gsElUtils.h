/** @file gsElUtils.h

    @brief Provides several simple utility and naming classes.

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

/// @brief Specifies the time integration scheme, see Wriggers, Nonlinear FEM, p. 205
struct time_integration
{
    enum scheme
    {
        explicit_ = 0,  /// explicit scheme
        explicit_lumped = 1,  /// explicit scheme with lumped mass matrix
        implicit_linear = 2,  /// implicit scheme with linear problem
        implicit_nonlinear = 3 /// implicit scheme with nonlinear problem
    };
};


/// @brief Specifies the verbosity of the Newton's solver
struct newton_verbosity
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2  /// output everything
    };
};

/// @brief Specifies the status of the Newton's solver
enum class newton_status { converged, /// method successfully converged
                           interrupted, /// solver was interrupted after exceeding the limit of iterations
                           working, /// solver working
                           bad_solution }; /// method was interrupted because the current solution is invalid

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

/// @brief Specifies the elasticity formulation: pure displacement or mixed displacement-pressure
enum class elasticity_formulation { displacement, mixed_pressure };

/** @brief Simple progress bar class
 *
 *         Prints the progress bar on a single console line, avoids clattering the console window and looks cool.
 *         Useful for programms with duration known in advance, e.g. transient simulations with a fixed number of time steps.
 *         Other console output will mess up the progress bar.
*/
class gsProgressBar
{
public:
    /// Constructor. Width is a number of symbols the progress bar spans
    gsProgressBar(index_t width = 25) : m_width(width) {}

    /// display the progress from 0 to 1
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
