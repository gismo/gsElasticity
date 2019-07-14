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

#include <gsElasticity/gsElBaseAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>

#include <gsElasticity/gsElUtils.h>

namespace gismo
{


// ToDo: -add second Piola-Kirchhoff stresses for nonlinear elasticity, both NeoHook and St.V.-K.
//       -add Neumann BC on the deformed configuration (currently Neumann BC is assumed to be set
//        in the reference configuration, dead-load problem)

/** @brief Assembles the stiffness matrix and the right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. The matrix and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
           Supports mixed displacement-pressure formulation.
*/
template <class T>
class gsElasticityAssembler : public gsElBaseAssembler<T>
{
public:
    typedef gsElBaseAssembler<T> Base;

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

    /// @brief Assembles the stiffness matrix and the RHS
    /// set *assembleMatrix* to false to only assemble the RHS;
    virtual void assemble(bool assembleMatrix = true);

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual bool assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix = true);

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for displacement formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual void assemble(const gsMultiPatch<T> & displacement, bool assembleMatrix = true);

    /// @ brief Assembles the tangential matrix and the residual for a iteration of Newton's method for mixed formulation;
    /// set *assembleMatrix* to false to only assemble the residual;
    /// ATTENTION: rhs() returns a negative residual (-r) !!!
    virtual void assemble(const gsMultiPatch<T> & displacement, const gsMultiPatch<T> & pressure,
                          bool assembleMatrix = true);

    /// @brief Construct displacement from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result) const;

    /// @brief Construct displacement and pressure from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T> & displacement, gsMultiPatch<T> & pressure) const;

    virtual void constructPressure(const gsMatrix<T> & solVector, gsMultiPatch<T> & pressure) const;


    /// @brief Construct Cauchy stress tensor for visualization (only valid for linear elasticity)
    void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                 gsPiecewiseFunction<T> & result,
                                 stress_type::type type = stress_type::von_mises) const;

    /** @brief Deform geometric domain using computed displacement field.
     *
     * Domain must have the same level of refinement as the basis used for computing the displacement.
     * Otherwise the number of control points will not match.
     */
    virtual void deformGeometry(const gsMatrix<T> & solVector, gsMultiPatch<T> & domain) const;

    /** @brief Set Dirichet degrees of freedom on a given side of a given patch from a given matrix.
     *
     * A usual source of degrees of freedom is another geometry where it is known that the corresponding
     * bases match. The function is not safe in that it assumes a correct numbering of DoFs in a given
     * matrix. To allocate space for these DoFs in the assembler, add an empty/zero Dirichlet boundary condition
     * to gsBoundaryCondtions container that is passed to the assembler constructor.
     */
    virtual void setDirichletDofs(size_t patch, boxSide side, const gsMatrix<T> & ddofs);

    /// @brief Check whether the displacement field is valid, i.e. J = det(F) > 0;
    /// return -1 if yes or a number of the first invalid patch
    virtual index_t checkSolution(const gsMultiPatch<T> & solution) const;

    /// @brief Return minJ/maxJ
    virtual T solutionJacRatio(const gsMultiPatch<T> & solution) const;

    /// sets scaling of Dirichlet BC used for linear system assembly
    virtual void setDirichletAssemblyScaling(T factor);
    /// sets scaling of Dirichlet BC used for construction of the solution as a gsMultiPatch object
    virtual void setDirichletConstructionScaling(T factor);
    /// set scaling of the force loading (volume and surface loading)
    virtual void setForceScaling(T factor);

protected:
    /// scale Dirichlet degrees of freedom
    void scaleDDoFs(T factor);
    /// reset Dirichlet degrees of freedom to its original state
    void resetDDoFs();

protected:

    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
	short_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
    /// Dirichlet degrees of freedom saved to recover after modification
    std::vector<gsMatrix<T> > saved_ddof;
};

/// @brief Generates a matrix of sampling points for a given parametric element;
/// includes quadrature points for the element as well as the corner points
template <class T>
void genSamplingPoints(const gsVector<T> & lower, const gsVector<T> & upper,
                       const gsQuadRule<T> & quRule, gsMatrix<T> & points);

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
