/** @file gsElasticityAssembler.h

    @brief Provides nonlinear elasticity system matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger, A.Shamanskiy (TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsElasticity/gsMultiFunction.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

namespace gismo
{

template< class T>
class gsStressFunction;

template< class T>
class gsMultiFunction;

struct stress_type { enum type{normal = 0, shear = 1, von_mises = 2}; };

/** @brief Assembles linear and non linear elasticity matrices for 2D plain strain and 3D continua.


    \tparam T coefficient type

    \ingroup Elasticity
*/
template <class T>
class gsElasticityAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

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
                            const gsFunction<T> & body_force,
                            dirichlet::values computeStrategy = dirichlet::l2Projection,
                            dirichlet::strategy enforceStrategy = dirichlet::elimination,
                            gsFunction<T> * E = nullptr,
                            gsFunction<T> * pr = nullptr);


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

    /** \brief Computes stresses \f$ \sigma_{ij}\f$ of
    * already computed solution.
    *
    * The computed stresses are written in a
    * Matrix of size <em>k</em> x <em>n</em>,
    * where <em>k=3</em> in 2D and <em>k=6</em> in 3D.\n
    *
    * If the von Mises stresses \f$ \sigma_v \f$, given by\n
    * in 2D:
    * \f[
      \sigma_v^2 =
      \sigma_{11}^2 + \sigma_{11}\sigma_{22} + \sigma_{22}^2 +
      3 \sigma_{12}^2,
    \f]
    * in 3D:
    * \f[
      \sigma_v^2 = \frac{1}{2}(
      (\sigma_{11}-\sigma_{22})^2 +
      (\sigma_{11}-\sigma_{33})^2 +
      (\sigma_{22}-\sigma_{33})^2 +
      6 ( \sigma_{12}^2 + \sigma_{13}^2 + \sigma_{23}^2 ) ),
    \f]
    * are also computed, these
    * are added as an additional line in the Matrix.
    *
    * Each column corresponds to an evaluation point and contains\n
    * in 2D, \em without von Mises stress:\n
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{12})^T\f$.\n
    * in 2D, \em with von Mises stress:\n
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{12},\sigma_v)^T\f$.\n
    * in 3D, \em without von Mises stress:\n
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{12},\sigma_{13},\sigma_{23})^T\f$.\n
    * in 3D, \em with von Mises stress:\n
    * \f$ (\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{12},\sigma_{13},\sigma_{23},\sigma_v)^T\f$.\n
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
    * computed stresses. See above for the format.
    * \param[in] computeVonMises indicates whether the von Mises stress
    * should also be computed. If \em true, the von Mises stress
    * is
    */
    void computeStresses(const gsMatrix<T>& solVector,
                         const gsMatrix<T>& u,
                         int patchIndex,
                         gsMatrix<T> &result,
                         bool computeVonMises = false) const;

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

    /// Generate patchwise stress-functions to be used for visualization.
    void constructStresses(const gsMatrix<T>& solVector,
                           gsMultiFunction<T>& result,
                           stress_type::type type) const;

    /// Add already computed neumann boundary data to the rhs vector.
    /// Resulting rhs is saved in rhsExtra member and can be accessed
    /// or cleared by the corresponding class methods.
    /// (orientation and size check for 2D/2.5D only)

    void addNeummannData(const gsMultiPatch<> & sourceGeometry,
                          const gsMultiPatch<> & sourceSolution,
                          int sourcePatch, const boxSide & sourceSide,
                          int targetPatch, const boxSide & targetSide);

    void addNeummannData(const gsField<> & sourceField,
                          int sourcePatch, const boxSide & sourceSide,
                          int targetPatch, const boxSide & targetSide);

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

    void clampPatchCorner(int patch,int corner);

    void deformGeometry(const gsMatrix<T> & solVector,
                        gsMultiPatch<T> & result);


    friend class gsStressFunction<T>;

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

    /// Check if two boundaries have at least the same starting and ending points
    /// Returns 1 if they match, return -1 if they match,
    /// but parameterization directions are opposite,
    /// return 0 if ends don't match
    int checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                const gsGeometry<> & targetBoundary);




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

    gsMatrix<T> m_rhsExtra;
    std::map<unsigned,T> externalDDofs;

    // Determines how the (fixed) Dirichlet values should be computed
    //dirichlet::values  m_dirValues;

    // Strategy for enforcing Dirichlet DoFs
    //dirichlet::strategy m_dirStrategy;

    gsFunction<T> * f_E;
    gsFunction<T> * f_pr;

protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

/** @brief Allows computation and visualization of von Mises stresses for linear elasticity.

    \tparam T coefficient type

    \ingroup Elasticity
*/
template <class T>
class gsStressFunction : public gsFunction<T>
{
public:

    gsStressFunction(const gsMatrix<T> & solVector, const gsElasticityAssembler<T> & assembler, index_t patchNum, stress_type::type type)
        : m_displacement(solVector),
          m_assembler(assembler),
          m_patch(patchNum),
          m_type(type)
    {

    }

    int domainDim() const
    {
        return m_assembler.m_bases[0].dim();
    }

    int targetDim() const
    {
        switch(m_type)
        {
        case stress_type::normal:
        {
            return domainDim();
        }
        case stress_type::shear:
        {
            return (domainDim() == 2) ? 1 : 3;
        }
        case stress_type::von_mises:
        {
            return 1;
        }
        default:
            gsInfo << "Bad stress type\n";
            return -1;
        };
    }

    void eval_into(const gsMatrix< T > & u,gsMatrix< T > & result ) const;


protected:
    gsMatrix<T> m_displacement;
    const gsElasticityAssembler<T>& m_assembler;
    index_t m_patch;
    stress_type::type m_type;

}; // class definition ends

} // namespace gismo



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityAssembler.hpp)
#endif
