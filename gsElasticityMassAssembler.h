/** @file gsElasticityMassAssembler.h

    @brief Provides elasticity system mass matrices for 2D plain strain and 3D continua.

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

/** @brief Assembles elasticity mass matrices for 2D plain strain and 3D continua.


    \tparam T coefficient type

    \ingroup Elasticity   
*/
template <class T>
class gsElasticityMassAssembler : public gsAssemblerBase<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    
    \ingroup Assembler
*/
    gsElasticityMassAssembler(  gsMultiPatch<T> const & patches,
		                    gsMultiBasis<T> const & bases,
							// material properties
							T density_rho,
							// Boundary conditions
							gsBoundaryConditions<T> const & bconditions
						);

public:

    /// Main assembly routine.
    void assemble();
    
protected:

    /// Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofsIntpl();

protected:

	T m_rho;

	/// Dimension (parameter space = physical space = deformation vector)
	index_t m_dim;

    /// Boundary conditions
    gsBoundaryConditions<T> m_bConditions;
 
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
    using gsAssemblerBase<T>::m_dofs;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityMassAssembler.hpp)
#endif
