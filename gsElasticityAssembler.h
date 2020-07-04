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

#include <gsElasticity/gsBaseAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>
#include <gsElasticity/gsBaseUtils.h>

namespace gismo
{

/** @brief Assembles the stiffness matrix and the right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. The matrix and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
           Supports mixed displacement-pressure formulation.
*/
template <class T>
class gsElasticityAssembler : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;

    /// @brief Constructor for displacement formulation
    gsElasticityAssembler(const gsMultiPatch<T> & patches,
                          const gsMultiBasis<T> & basis,
                          const gsBoundaryConditions<T> & bconditions,
                          const gsFunction<T> & body_force);

    /// @brief Constructor of mixed formulation (displacement + pressure)
    gsElasticityAssembler(const gsMultiPatch<T> & patches,
                          const gsMultiBasis<T> & basisDisp,
                          const gsMultiBasis<T> & basisPres,
                          const gsBoundaryConditions<T> & bconditions,
                          const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    //--------------------- SYSTEM ASSEMBLY ----------------------------------//

    /// @brief Assembles the stiffness matrix and the RHS for the LINEAR ELASTICITY
    /// set *assembleMatrix* to false to only assemble the RHS;
    virtual void assemble(bool saveEliminationMatrix = false);

    /// Assembles the tangential linear system for Newton's method given the current solution
    /// in the form of free and fixed/Dirichelt degrees of freedom.
    /// Checks if the current solution is valid (Newton's solver can exit safely if invalid).
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs);
protected:
    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for displacement formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual void assemble(const gsMultiPatch<T> & displacement);

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for mixed formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual void assemble(const gsMultiPatch<T> & displacement, const gsMultiPatch<T> & pressure);

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//

public:
    /// @brief Construct displacement from computed solution vector and fixed degrees of freedom
    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & displacement) const;

    /// @brief Construct displacement and pressure from computed solution vector and fixed degrees of freedom
    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & displacement, gsMultiPatch<T> & pressure) const;

    /// @ brief Construct pressure from computed solution vector
    virtual void constructPressure(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & pressure) const;

    //--------------------- SPECIALS ----------------------------------//

    /// @brief Construct Cauchy stresses for evaluation or visualization
    virtual void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                 gsPiecewiseFunction<T> & result,
                                 stress_components::components component = stress_components::von_mises) const;

    /// @brief Construct Cauchy stresses for evaluation or visualization
    virtual void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                 const gsMultiPatch<T> & pressure,
                                 gsPiecewiseFunction<T> & result,
                                 stress_components::components component = stress_components::von_mises) const;

protected:
    /// a custom reserve function to allocate memory for the sparse matrix
    virtual void reserve();

protected:
    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
    short_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
    using Base::eliminationMatrix;
};


} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
