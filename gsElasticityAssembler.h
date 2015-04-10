/** @file gsElasticityAssembler.h

    @brief Provides nonlinear elasticity system matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once

#include <gsAssembler/gsAssemblerBase.h>

namespace gismo
{

/** @brief Assembles linear and non linear elasticity matrices for 2D plain strain and 3D continua.


    \tparam T coefficient type

    \ingroup Elasticity   
*/
template <class T>
class gsElasticityAssembler : public gsAssemblerBase<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    
    \ingroup Assembler
*/
    gsElasticityAssembler(  gsMultiPatch<T> const & patches,
		                    gsMultiBasis<T> const & bases,
							// material properties
							T E_modulus,
							T poissons_ratio,
							T density_rho,
							// Boundary conditions
							gsBoundaryConditions<T> const & bconditions,
							// Body force per unit surface/volume (in 2D/3D)
							const gsFunction<T> & body_force
							// Points on physical domain
							// const gsMatrix<T> pointCoords,
							// Point forces
							// const gsMatrix<T> pointForces
							);

public:

    /// Main assembly routine.
    void assemble();

    /// Main assembly routine for the non-linear case
    void assemble(const gsMultiPatch<T> & deformed);

    /// Reconstruct solution from computed solution vector
    void constructSolution(const gsMatrix<T>& solVector, 
                           gsMultiPatch<T>& result) const;

    // Set solution from solVector, overwrites previous solution
    void setSolution(const gsMatrix<T>& solVector, 
                     gsMultiPatch<T>& result) const;

	 // Newton update of the solution from solVector
    void updateSolution(const gsMatrix<T>& solVector, 
                        gsMultiPatch<T>& result) const;

	/// Set factor for time-dependent external forces (at current time-step)
    void set_tfac(const T tfac_neumann,
		          const T tfac_force);
    
protected:

    /// Neumann contributions
    void assembleNeumann();

    /// Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofsIntpl();

protected:

	/// Material parameters
    T m_lambda;
    T m_mu;
	T m_rho;

	/// Dimension (parameter space = physical space = deformation vector)
	index_t m_dim;

    /// Boundary conditions
    gsBoundaryConditions<T> m_bConditions;

    const gsFunction<T> *m_bodyForce;

	/// Factor for time-dependent external forces
	T m_tfac_neumann;
	T m_tfac_force;
   
    // Determines how the (fixed) Dirichlet values should be computed
    //dirichlet::values  m_dirValues;

    // Strategy for enforcing Dirichlet DoFs
    //dirichlet::strategy m_dirStrategy;

protected:

    // Members from gsAssemblerBase
    using gsAssemblerBase<T>::m_patches;
    using gsAssemblerBase<T>::m_bases;
    using gsAssemblerBase<T>::m_dofMappers;
    using gsAssemblerBase<T>::m_ddof;
    using gsAssemblerBase<T>::m_matrix;
    using gsAssemblerBase<T>::m_rhs;
    using gsAssemblerBase<T>::m_dofs;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
