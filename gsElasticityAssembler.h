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

#include <gsAssembler/gsAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

namespace gismo
{

/** @brief Assembles stiffness and mass matrices and right-hand side vector for linear and nonlinear elasticity
           for 2D plain stress and 3D continua. Matrices and vector have a block structure associated with
           components of the displacement vector, each block corresponding to one component.
*/
template <class T>
class gsElasticityAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

    /// @brief Constructor of the assembler object.
    gsElasticityAssembler(  gsMultiPatch<T> const & patches,
                            gsMultiBasis<T> const & bases,
                            gsBoundaryConditions<T> const & bconditions,
                            const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the stiffness matrix
    virtual void assemble();

    /// @brief Construct patchwise stress-function for visualization.
    void constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                 gsPiecewiseFunction<T>& result,
                                 stress_type::type type = stress_type::von_mises) const;

    //=========================================//
    //============ Old/to-do functions ========//
    //=========================================//

    /// @brief Changes the coefficient of an already existing multipatch
    /// so that they correspond to a given solution vector
    void setSolution(const gsMatrix<T>& solVector,
                     gsMultiPatch<T>& result) const;

    /// @brief Sum the coefficient of an already existing multipatch with
    /// those corresponding to a given solution vector
    void updateSolution(const gsMatrix<T>& solVector,
                        gsMultiPatch<T>& result) const;

    /// @brief (NonLin) Recompute Dirichlet DoFs for a deformed configuration
	void reComputeDirichletDofs(gsMultiPatch<T> &deformed);


    /// Add already computed Dirichlet boundary data to the specified sides.
    /// Orientation check only for 2D.
    void addDirichletData(const gsMultiPatch<> & sourceGeometry,
                          const gsMultiPatch<> & sourceSolution,
                          int sourcePatch, const boxSide & sourceSide,
                          int targetPatch, const boxSide & targetSide);

    void addDirichletData(const gsField<> & sourceField,
                          int sourcePatch, const boxSide & sourceSide,
                          int targetPatch, const boxSide & targetSide);

    void setDirichletDoFs(const gsMatrix<> & ddofs,
                          int targetPatch,
                          const boxSide & targetSide);


    const gsMatrix<T> & rhs();

    void resetRhsExtra() { m_rhsExtra.clear(); }

    void deformGeometry(const gsMatrix<T> & solVector,
                        gsMultiPatch<T> & result);

    /// @brief (NonLin) Assembles the stiffness matrix for a given deformed configuration
    void assemble(const gsMultiPatch<T> & deformed);

    /// @brief Assembles the mass matrix for the eigen-value analysis
    void assembleMass();


protected:

    int checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                const gsGeometry<> & targetBoundary);

protected:

    /// Dimension of the problem
    /// parametric dim = physical dim = deformation dim
	index_t m_dim;

    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

    //=========================================//
    //============ Old/to-do members ==========//
    //=========================================//

    gsMatrix<T> m_rhsExtra;
    std::map<unsigned,T> externalDDofs;
};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
