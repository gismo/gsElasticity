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

    /// @brief (NonLin) Assembles the stiffness matrix for a given deformed configuration
    void assemble(const gsMultiPatch<T> & deformed);

    /// @brief Assembles the mass matrix for the eigen-value analysis
	void assembleMass();

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

    /// @brief Construct patchwise stress-functions for visualization.
    void constructStresses(const gsMatrix<T>& solVector,
                           gsMultiFunction<T>& result,
                           stress_type::type type) const;


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

    friend class gsStressFunction<T>;

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

    gsMatrix<T> m_rhsExtra;
    std::map<unsigned,T> externalDDofs;
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
