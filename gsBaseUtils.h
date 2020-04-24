/** @file gsBaseUtils.h

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
#include <gsUtils/gsUtils.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{

/// @brief Specifies method used for mesh deformation in fluid-structure interaction
struct ale_method
{
    enum method
    {
        HE = 0,         /// harmonic extension
        IHE = 1,        /// incremental harmonic extension
        LE = 2,         /// linear elasticity
        ILE = 3,        /// incremental linear elasticity
        TINE = 4,       /// tangential incremental nonlinear elasticity (with the neo-Hookean law)
        TINE_StVK = 5,  /// tangential incremental nonlinear elasticity (with the St.Venant-Kirchhoff law)
        BHE = 6,        /// bi-harmonic extension
        IBHE = 7,       /// incremental bi-harmonic extension
    };
};

/// @brief Specifies the iteration type used to solve nonlinear systems
struct ns_assembly
{
    enum type
    {
        ossen = 0,          /// stationary point iteration, 1st order, yields a new solution at each iteration
        newton_update = 1,  /// newton's method, 2nd order, yields updates to the solution
        newton_next = 2     /// newton's method, 2nd order, yields a new solution at each iteration
    };
};

/// @brief Specifies the time integration scheme, see Wriggers, Nonlinear FEM, p. 205
struct time_integration
{
    enum scheme
    {
        explicit_ = 0,         /// explicit scheme
        explicit_lumped = 1,   /// explicit scheme with lumped mass matrix
        implicit_linear = 2,   /// implicit scheme with linear problem (theta-scheme)
        implicit_nonlinear = 3 /// implicit scheme with nonlinear problem (theta-scheme)
    };
};

/// @brief Specifies linear solver to use if it is hidden within some other class (like Newton's method or time integrators)
struct linear_solver
{
    enum solver
    {
        LU = 0,              /// LU decomposition: direct, no matrix requirements, robust but a bit slow, Eigen and Pardiso available
        LDLT = 1,            /// Cholesky decomposition pivoting: direct, simmetric positive or negative semidefinite, rather fast, Eigen and Pardiso available
        CGDiagonal = 2,      /// Conjugate gradient solver with diagonal (a.k.a. Jacobi) preconditioning: iterative(!), simmetric, Eigen only
        BiCGSTABDiagonal = 3 /// Bi-conjugate gradient stabilized solver with diagonal (a.k.a. Jacobi) preconditioning: iterative(!), no matrix requirements, Eigen only
    };
};


/// @brief Specifies the verbosity of the iterative solver
struct solver_verbosity
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2    /// output everything
    };
};

/// @briefly Specifies iteration type for an iterative solver
struct iteration_type
{
    enum type
    {
        update = 0, /// each iteration yields an update
        next = 1    /// each iteration yields a next solution
    };
};

/// @brief Specifies the status of the iterative solver
enum class solver_status { converged,      /// method successfully converged
                           interrupted,    /// solver was interrupted after exceeding the limit of iterations
                           working,        /// solver working
                           bad_solution }; /// method was interrupted because the current solution is invalid

/** @brief Specifies the type of stresses to compute
 *
 *         Currently, gsWriteParaview can only plot vector-valued functions with an output dimension up to three.
 *         Therefore it not possible to plot all stress components as components of a single vector-valued function.
*/
struct stress_components
{
    enum components
    {
        von_mises = 0,         /// return von Mises stress as a scala
        all_2D_vector    = 1,  /// return all components of the 2D stress tensor as a 3x1 vector
                               /// (s11 s22 s12) (useful for Paraview plotting)
        all_2D_matrix = 2,     /// return all components of the 2D stress tensor as a 2x2 matrix
        normal_3D_vector = 3,  /// return normal components of the 3D stress tensor as a 3x1 vector
                               /// (s11 s22 s33) (useful for Paraview plotting)
        shear_3D_vector  = 4,  /// return shear components of the 3D stress tensor as a 3x1 vector
                               /// (s12 s13 s23) (useful for Paraview plotting)
        all_3D_matrix = 5      /// return all components of the 3D stress tensor as a 3x3 matrix
    };
};

/// @brief Specifies the material law to use
struct material_law
{
    enum law
    {
        mixed_hooke               = -2, /// sigma = 2*mu*eps + p*I
        hooke                     = -1, /// sigma = 2*mu*eps + lambda*tr(eps)*I
        saint_venant_kirchhoff    = 0, /// S = 2*mu*E + lambda*tr(E)*I
        neo_hooke_ln              = 1, /// S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
        neo_hooke_quad            = 2, /// S = lambda/2*(J^2-1)C^-1 + mu*(I-C^-1)
        mixed_neo_hooke_ln        = 3 /// S = p*C^-1 + mu*(I-C^-1)
        //mixed_kelvin_voigt        = 4 /// S = p*C^-1 + mu*(I-C^-1*tr(C)/3)+ nu*C^-1*C'*C^-1
        //mixed_kelvin_voigt_muscle = 5  /// S = p*C^-1 + gamma*[ mu_m*(I-C^-1*tr(C)/3)+ nu_m*C^-1*C'*C^-1 ]
        //                                      + (1-gamma)*[ mu_t*(I-C^-1*tr(C)/3)+ nu_t*C^-1*C'*C^-1 ]
        // here, gamma in [0,1] is an indicator of muscle tissue; (1-gamma) indicates tendon tissue
    };
};

struct GISMO_EXPORT gsBoundaryInterface
{
    gsBoundaryInterface() {}

    /// boundary interface sides
    std::vector<patchSide> sidesA;
    std::vector<patchSide> sidesB;
    /// patch-to-patch correspondence
    std::vector<std::pair<index_t,index_t> > patches;


    void addInterfaceSide(index_t patchA, boundary::side sideA,
                          index_t patchB, boundary::side sideB)
    {
        sidesA.push_back(patchSide(patchA,sideA));
        sidesB.push_back(patchSide(patchB,sideB));
    }

    void addPatches(index_t patchA, index_t patchB)
    { patches.push_back(std::pair<index_t,index_t>(patchA,patchB)); }
};

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
        index_t threshold = index_t(progress*m_width);
        gsInfo << "[";
        for(index_t i = 0; i < m_width; i++)
            if(i < threshold)
                gsInfo << "=";
            else if(i == threshold)
                gsInfo << ">";
            else
                gsInfo << " ";
        gsInfo << "] " << (abs(progress - 1.) < 1e-12 ? 100 : index_t(progress*100)) << " %\r";
        gsInfo.flush();

        if (abs(progress - 1.) < 1e-12)
            gsInfo << std::endl;
    }

    /// display the progress from 0 to 1
    void display(index_t progress, index_t total)
    {
        GISMO_ENSURE(progress >= 0 && progress <= total && total >= 0,"Invalid progress value!");
        index_t threshold = index_t(1.*progress*m_width/total);
        gsInfo << "[";
        for(index_t i = 0; i < m_width; i++)
            if(i < threshold)
                gsInfo << "=";
            else if(i == threshold)
                gsInfo << ">";
            else
                gsInfo << " ";
        gsInfo << "] " << progress << "/" << total << " \r";
        gsInfo.flush();

        if (progress == total)
            gsInfo << std::endl;
    }

protected:
    index_t m_width;
};

template <class T>
std::string secToHMS(T sec)
{
    if (sec < 10)
        return util::to_string(sec) + "s";

    index_t days = index_t(sec)/(3600*24);
    index_t residual = index_t(sec)- 3600*24*days;
    index_t hours = residual/3600;
    residual -= 3600*hours;
    index_t minutes = residual/60;
    residual -= 60*minutes;

    std::string result = util::to_string(residual) + "s";
    if (minutes > 0)
        result = util::to_string(minutes) + "m" + result;
    if (hours > 0)
        result = util::to_string(hours) + "h" + result;
    if (days > 0)
        result = util::to_string(days) + "d" + result;

    return result;
}

}
