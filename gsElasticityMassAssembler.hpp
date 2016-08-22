/** @file gsElasticityMassAssembler.hpp

    @brief Provides elasticity system mass matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once

#include <gsElasticity/gsElasticityMassAssembler.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorMassElasticity.h>

// ---
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo
{

template<class T>
gsElasticityMassAssembler<T>::gsElasticityMassAssembler(gsMultiPatch<T> const & patches,
														gsMultiBasis<T> const & bases,
														T density_rho,
														gsBoundaryConditions<T> const & bconditions
														)
:   Base(patches), 
	m_rho(density_rho),
    m_bConditions(bconditions)
{
    // Copy the spline basis
    m_bases.push_back(bases);
       
	m_dim = m_patches.parDim();
	GISMO_ASSERT(m_dim == m_patches.geoDim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == m_bases.front().dim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only for 2D and 3D!");
	
    // Mark boundary dofs
    m_dofMappers.resize(m_dim);
	for (index_t i = 0; i < m_dim; i++)
		m_bases.front().getMapper(true, m_bConditions, m_dofMappers[i] ); 

    // Determine system size
    m_dofs = 0;
	for (index_t i = 0; i < m_dim; i++)
        m_dofs += m_dofMappers[i].freeSize();

    // Add shifts to global matrix indices to get the global dof
    // number for each unknown coordinate
	unsigned inShift = 0;
	unsigned bdShift = 0;

	for (index_t i = 0; i < m_dim; i++)
	{
        m_dofMappers[i].setShift( inShift );
		m_dofMappers[i].setBoundaryShift( bdShift );
		inShift += m_dofMappers[i].freeSize();
		bdShift += m_dofMappers[i].boundarySize();
	}
}


template<class T>
void gsElasticityMassAssembler<T>::assemble()
{
    std::cout << "Linear Elasticity: assemble mass matrix." << std::endl;
	
	//computeDirichletDofs();
    index_t numDirichlet = 0;
    for (index_t i = 0; i < m_dim; ++i)
        numDirichlet += m_dofMappers[i].boundarySize();
    m_ddof.setZero(numDirichlet, 1);

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    size_t nonZerosPerCol = m_dim;
    for (index_t i = 0; i < m_dim; ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
        
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1 );

    // Assemble volume stiffness and load vector integrals
    gsVisitorMassElasticity<T> visitor(m_rho);
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T> // AM: is non-homogeneous not called yet
void gsElasticityMassAssembler<T>::computeDirichletDofsIntpl()
{
    const gsDofMapper &mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), 1 ); //--mrhs
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); 
		  ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary->size(); ++i)
            {
                const int ii= mapper.bindex( (*boundary)(i) , k );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( m_dim );

        for ( int i = 0; i < m_dim; ++i)
        {
            if ( i == dir )
            {
                gsVector<T> b(1); 
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( basis.component(i).anchors()->transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<T> fpts = 
            it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid<T>( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        gsGeometry<T> * geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k2=0; k2!= boundary->size(); ++k2)
        {
            const int ii= mapper.bindex( (*boundary)(k2) , it->patch() );
            m_ddof.row(ii) = dVals.row(k2);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}


}// namespace gismo
