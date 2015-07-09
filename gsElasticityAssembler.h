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

    /** \brief Main assembly routine.
     *
     * Assembles the system matrix. The vector of unknowns
     * is sorted such that all DOFs associated with the
     * first coordinate direction of the solution
     * come first, then all
     * DOFs associated with the second coordinate direction,
     * and last, all DOFs associated with the third
     * coordinate direction.
     */
    void assemble();

    /// Main assembly routine for the non-linear case
    void assemble(const gsMultiPatch<T> & deformed);

	void assembleMass();

    /** \brief Reconstruct solution from computed solution vector
     *
     * \param[in] solVector The vector of coefficients for the
     * \em free DOF (i.e., those that are NOT on the
     * Dirichlet (displacement) boundary.\n
     * Given as
     * gsMatrix of size <em>d*n</em> x \em 1, where\n
     * \em d is the dimension of the domain and\n
     * \em n is the number of free DOF.\n
     * The coefficients are ordered such that (in 3D) the
     * first \em n coefficients correspond to basis
     * functions of type \f$ (B_i,0,0)^T\f$, the next \em n
     * coefficients corrspond to basis
     * functions of type \f$ (0,B_i,0)^T\f$, and so on.
     * \param[out] result
     */
    void constructSolution(const gsMatrix<T>& solVector, 
                           gsMultiPatch<T>& result) const;

	/** \brief Computes stresses \f$ \sigma_{ij}\f$ of    * already computed solution.
    *
    * \param[in] solVector Solution vector containing the
    * computed \em free degrees of freedom (i.e., without
    * the DOFs on the Dirichlet (=displacement) boundary
    * \param[in] u Evaluation points as matrix of size
    * <em>d</em> x <em>n</em> where \em d is the dimension
    * of the domain and \em n the number of evaluation points.
    * Each column in \em u corresponds to one point in the
    * parameter domain.
    * \param[in] patchIndex Index of the patch on which
    * the computation shall be done.
    * \param[out] result Gets overwritten in with the
    * computed stresses. Matrix of size <em>k</em> x <em>n</em>,
    * where <em>k=3</em> in 2D and <em>k=6</em> in 3D.\n
    * Each column corresponds to an evaluation point and contains\n
    * in 2D:
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{12})^T\f$.\n
    * in 3D:
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{12},\sigma_{13},\sigma_{23})^T\f$.\n
    */
    void computeStresses(const gsMatrix<T>& solVector,
                         const gsMatrix<T>& u,
                         int patchIndex,
                         gsMatrix<T> &result) const;

    // Set solution from solVector, overwrites previous solution
    void setSolution(const gsMatrix<T>& solVector, 
                     gsMultiPatch<T>& result) const;

	 // Newton update of the solution from solVector
    void updateSolution(const gsMatrix<T>& solVector, 
                        gsMultiPatch<T>& result) const;

	/// Set factor for time-dependent external forces (at current time-step)
    void set_tfac(const T tfac_neumann,
		          const T tfac_force);

	/// Set the \em m_MATERIAL_LAW for nonlinear assembly to 0: St. Venant-Kirchhoff, 1: Neo-Hooke
    void set_MaterialLaw(const int material);

	/// Re-Compute Dirichlet DoFs after Update and set deformed to correct values
	///   needed for nonlinear with changing Dirichlet BC (displacement control)
	void reComputeDirichletDofs(gsMultiPatch<T> &deformed);
    
protected:

    /// Neumann contributions
    void assembleNeumann();

	/** Computes the Dirichlet DoF values by interpolation 
	    * 
		* \warning Works only for tensor-product-bases! 
	*/ 
    void computeDirichletDofsIntpl();

	/** \brief Computes Dirichlet-boundary conditions by L2-projection.
		*
		* ...if the dirichlet::strategy is chosen as dirichlet::elimination.\n
		* A global \f$L_2\f$-projection is applied to the given Dirichlet
		* (displacement) data
		* and the eliminated coefficients are set to the corresponding values.
		* The projection is global in the sense that all Dirichlet-DOFs are
		* computed at once.
		*
	*/
	void computeDirichletDofsL2Proj();

protected:

	/// Material parameters
    T m_lambda;
    T m_mu;
	T m_rho;
	int m_MATERIAL_LAW;

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
