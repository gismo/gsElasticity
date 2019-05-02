/** @file gsMixedElasticityAssembler.h

    @brief Provides linear and nonlinear mixed elasticity systems for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

namespace gismo
{

template <class T>
class gsMixedElasticityAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

    /// @brief Constructor of the assembler object.
    gsMixedElasticityAssembler(const gsMultiPatch<T> & patches,
                               const gsMultiBasis<T> & basisDisplacement,
                               const gsMultiBasis<T> & basisPressure,
                               const gsBoundaryConditions<T> & bconditions,
                               const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the stiffness matrix and the RHS
    virtual void assemble();

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    //virtual void assemble(const gsMultiPatch<T> & deformed, bool assembleMatrix = true);

    /// @brief Construct solution from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T> & displacement, gsMultiPatch<T> & pressure) const;

    /// @brief Construct Cauchy stress tensor for visualization (only valid for linear elasticity)
    //void constructCauchyStresses(const gsMultiPatch<T> & displacement,
    //                             gsPiecewiseFunction<T> & result,
    //                             stress_type::type type = stress_type::von_mises) const;

    /** @brief Deform geometric domain using computed displacement field.
     *
     * Domain must have the same level of refinement as the basis used for computing the displacement.
     * Otherwise the number of control points will not match.
     */
    //virtual void deformGeometry(const gsMatrix<T> & solVector, gsMultiPatch<T> & domain) const;

    /** @brief Set Dirichet degrees of freedom on a given side of a given patch from a given matrix.
     *
     * A usual source of degrees of freedom is another geometry where it is known that the corresponding
     * bases match. The function is not safe in that it assumes a correct numbering of DoFs in a given
     * matrix. To allocate space for these DoFs in the assembler, add an empty/zero Dirichlet boundary condition
     * to gsBoundaryCondtions container that is passed to the assembler constructor.
     */
    //virtual void setDirichletDofs(index_t patch, boxSide side, const gsMatrix<T> & ddofs);

    /// @brief Check whether the displacement field is valid, i.e. J = det(F) > 0;
    /// return -1 if yes or a number of the first invalid patch
    //virtual index_t checkSolution(const gsMultiPatch<T> & solution) const;

    /// @brief Return minJ/maxJ
    //virtual T solutionJacRatio(const gsMultiPatch<T> & solution) const;

protected:

    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
    short_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

};

/// @brief Generates a matrix of sampling points for a given parametric element;
/// includes quadrature points for the element as well as the corner points
//template <class T>
//void genSamplingPoints(const gsVector<T> & lower, const gsVector<T> & upper,
//                       const gsQuadRule<T> & quRule, gsMatrix<T> & points);

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
