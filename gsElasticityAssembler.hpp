/** @file gsElasticityAssembler.hpp

    @brief Provides nonlinear elasticity system matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger, A. Shamanskiy (TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
#include <gsElasticity/gsVisitorNonLinElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>
#include <gsElasticity/gsVisitorElasticityPressure.h>
#include <gsElasticity/gsVisitorMassElasticity.h>

// ---
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
												gsMultiBasis<T> const & bases,
												T E_modulus,
												T poissons_ratio,
												T density_rho,
												gsBoundaryConditions<T> const & bconditions,
                                                const gsFunction<T> & body_force,
                                                dirichlet::values computeStrategy,
                                                dirichlet::strategy enforceStrategy)
:   Base(patches), 
	m_rho(density_rho),
    m_bConditions(bconditions),
    m_bodyForce(&body_force)
{
    // Copy the spline basis
    m_bases.push_back(bases);
    
    // Initialize material properties
    m_lambda = E_modulus * poissons_ratio / ( (1.+poissons_ratio)*(1.-2.*poissons_ratio)) ;
    m_mu     = E_modulus / (2.*(1.+poissons_ratio)) ;
    
	m_dim = body_force.targetDim();
	GISMO_ASSERT(m_dim == m_patches.parDim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == m_bases.front().dim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only for 2D and 3D!");
	
    // Mark boundary dofs
    m_dofMappers.resize(m_dim);
	for (index_t i = 0; i < m_dim; i++)
		m_bases.front().getMapper(true, m_bConditions, i, m_dofMappers[i], true); // m_bases.front().getMapper(true, m_bConditions, m_dofMappers[i]); 

    // Determine system size
    m_dofs = 0;
	for (index_t i = 0; i < m_dim; i++)
        m_dofs += m_dofMappers[i].freeSize();

	// Initialize factors for time-dependent external forces
	m_tfac_neumann = 1.0;
	m_tfac_force = 1.0;

	m_MATERIAL_LAW = 1; 


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

    m_options.dirValues = computeStrategy;
    m_options.dirStrategy = enforceStrategy;

    computeDirichletDofs();
}

template<class T>
void gsElasticityAssembler<T>::set_tfac(const T tfac_neumann,
		                                const T tfac_force)
{
	m_tfac_neumann = tfac_neumann;
	m_tfac_force = tfac_force;
}

/// Sets the \em m_MATERIAL_LAW to 0: St. Venant-Kirchhoff, 1: Neo-Hooke
template<class T>
void gsElasticityAssembler<T>::set_MaterialLaw(const int material)
{
    GISMO_ASSERT( (material >= 0) && (material <= 1), "Only 0 (St. Venant-Kirchhoff) or 1 (Neo-Hooke) allowed!");
	m_MATERIAL_LAW = material;
}

template<class T>
void gsElasticityAssembler<T>::assembleNeumann()
{
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd(); 
		  ++it )
    {
        gsVisitorElasticityNeumann<T> neumann(*it->function(), it->side(), m_tfac_neumann);

        // Note: it->unknown()
        this->apply(neumann, it->patch(), it->side() );
    }
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{

    //computeDirichletDofs();

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
    m_rhs.resize(m_dofs, 1 );
    m_rhs.setZero();

    // Assemble volume stiffness and load vector integrals
    gsVisitorLinearElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, m_tfac_force);
    for (unsigned np = 0; np < m_patches.nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann boundary conditions
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
    std::cout << "Nonlinear elasticity: assemble residual and tangential matrix." << std::endl;
	
	if ( m_ddof.size() == 0 )
    {
        // assemble(); and uncomment "computeDirichletDofs();" in assemble(); and comment "computeDirichletDofs();" in constructor
        computeDirichletDofs();
        return;
    }


    // Initialize the matrix and rhs vector
    m_matrix.setZero();
    
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1);

    gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0), m_tfac_force);
	//gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0));

	visitor.set_MaterialLaw(m_MATERIAL_LAW);

    // Assemble volume stiffness and load vector integrals
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        visitor.setDeformed( deformed.patch(np) );

        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann forces
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

size_t sumBoundarySize(size_t a, gsDofMapper b) { return a + b.boundarySize(); }

template<class T>
void gsElasticityAssembler<T>::computeDirichletDofs()
{
    size_t numDDofs = std::accumulate(m_dofMappers.begin(), m_dofMappers.end(), (size_t)0, sumBoundarySize);
    m_ddof.setZero(numDDofs,1);

    if (m_options.dirValues == dirichlet::interpolation)
    {
        computeDirichletDofsIntpl();
    }
    else if (m_options.dirValues == dirichlet::l2Projection)
    {
        computeDirichletDofsL2Proj();
    }
    else
    {
        gsInfo << "This Dirichlet value strategy is not implemented\n";
    }
}


// WARNING: WORKS ONLY FOR TENSOR-PRODUCT-BASES!
template<class T> 
void gsElasticityAssembler<T>::computeDirichletDofsIntpl()
{   
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); 
		  ++it )
    {
        if ( it->isHomogeneous() )
			continue;

		const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk/m_dim])[k];

		const gsDofMapper &mapper = m_dofMappers[unk];

        // Get dofs on this boundary
        gsMatrix<unsigned> boundary = basis.boundary(it->side()) ;

		/*
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
		*/

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
                rr.push_back( basis.component(i).anchors().transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<T> fpts = 
            it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid<T>( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k2=0; k2!= boundary.size(); ++k2)
        {
            const int ii= mapper.bindex( (boundary)(k2) , it->patch() );
            m_ddof.row(ii) = m_tfac_neumann * dVals.row(k2);
        }
        delete h;
    }
}

template<class T>
void gsElasticityAssembler<T>::computeDirichletDofsL2Proj()
{
    size_t totalDirDofs = 0;
    for( size_t i=0; i < size_t( m_dofMappers.size() ); i++)
        totalDirDofs += m_dofMappers[i].boundarySize();
    // ...note that gsDofMapper.boundarySize() counts only to the
    // eliminated DOFs.

    m_ddof.resize( totalDirDofs, 1 );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs( totalDirDofs, 1 );
    globProjRhs.setZero();

    // Temporaries
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals(1,1);
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals;

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_bConditions.dirichletBegin();
          iter != m_bConditions.dirichletEnd(); ++iter )
    {
        if ( iter->isHomogeneous() )
            continue;

        const int unk = iter->unknown();
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[0])[patchIdx];

        typename gsGeometry<T>::Evaluator geoEval( m_patches[patchIdx].evaluator(NEED_MEASURE));

        // Set up quadrature. The number of integration points in the direction
        // that is NOT along the element has to be set to 1.
        gsVector<index_t> numQuNodes( basis.dim() );
        for( int i=0; i < int( basis.dim() ); i++)
            numQuNodes[i] = (basis.degree(i)+1);
        numQuNodes[ iter->side().direction()] = 1;

        gsGaussRule<T> bdQuRule(numQuNodes);

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                          quNodes, quWeights);

            geoEval->evaluateAt( quNodes );

            // the values of the boundary condition are stored
            // to m_rhsGradV. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            gsMatrix<T> quPhys = m_patches[patchIdx].eval( quNodes );
            rhsVals.resize( quPhys.rows(), quPhys.cols() );
            iter->function()->eval_into( quPhys , rhsVals );

            basis.eval_into( quNodes, basisVals);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(quNodes.col(0), globIdxAct );
            m_dofMappers[unk].localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(m_dofMappers[unk].boundarySize());
            for( size_t i=0; i < size_t( globIdxAct.rows() ); i++)
                if( m_dofMappers[unk].is_boundary_index( globIdxAct(i,0)) )
                {
                    eltBdryFcts.push_back( i );
                    //isB( globIdxAct(i,0) , unk ) = globIdxAct(i,0);
                }

            // Do the actual assembly:
            for( index_t k=0; k < quNodes.cols(); k++ )
            {
                const T weight_k = quWeights[k] * geoEval->measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < size_t( eltBdryFcts.size() ); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = m_dofMappers[unk].global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < size_t( eltBdryFcts.size() ); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = m_dofMappers[unk].global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( totalDirDofs, totalDirDofs );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_ddof = solver.compute( globProjMat ).solve ( globProjRhs );

} // computeDirichletDofsL2Proj


template<class T>
void  gsElasticityAssembler<T>::reComputeDirichletDofs(gsMultiPatch<T> &deformed)
{
	computeDirichletDofsIntpl();

	// -> deformed must be set using updated m_ddof! How?
	for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update displacement solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = deformed.patch(p).coefs();
	
        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
				if ( !mapper.is_free(i, p) ) // eliminated DoF: fill with Dirichlet data              
                {
                    coeffs(i,j) = m_ddof(mapper.bindex(i, p), 0);
                } 
            }
        }
	}

}


template<class T>
void  gsElasticityAssembler<T>::setSolution(const gsMatrix<T>& solVector, 
                                            gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();
		coeffs.setZero(); //OWEE

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
				if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof(mapper.bindex(i, p), 0);
                }
            }
        }
    }
}


template<class T>
void  gsElasticityAssembler<T>::updateSolution(const gsMatrix<T>& solVector, 
                                               gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
            }
        }
    }
}

template<class T>
void  gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector, 
                                                  gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");
    
    // empty the multipatch container
    result.clear();

	gsMatrix<T> coeffs;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size(); // m_patches[p].size();

        coeffs.resize(sz, m_dim);

        for (index_t j = 0; j < m_dim; ++j) // For all components x, y, z
        {
            const gsDofMapper & mapper = m_dofMappers[j];// grab mapper for this component

            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) = solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof(mapper.bindex(i, p), 0);
                }
            }
        }

		result.addPatch( m_bases[0][p].makeGeometry( give(coeffs) ) );
    }
}


template<class T>
void gsElasticityAssembler<T>::assembleMass()
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


template<class T>
void gsElasticityAssembler<T>::computeStresses(
        const gsMatrix<T>& solVector,
        const gsMatrix<T>& u,
        int patchIndex,
        gsMatrix<T>& result,
        bool computeVonMises ) const
{
    //m_dim = m_basis.dim();
    size_t dimStrain = (m_dim*(m_dim+1))/2;

    result.resize( dimStrain + unsigned( computeVonMises), u.cols() );
    result.setZero();

    unsigned evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    typename gsGeometry<T>::Evaluator geoEval(
        m_patches[patchIndex].evaluator(evFlags));

    // Evaluate basis functions on element
    gsMatrix<T> basisGrad( m_dim, 1 );
    gsMatrix<unsigned> actives;

    gsMatrix<T> physGrad( m_dim, 1 );
    gsMatrix<T> uPartDers( m_dim, m_dim );

    // Do this pointwise, since we have no assumptions
    // on the evaluation points.
    for( size_t k=0; k < size_t( u.cols() ); k++ )
    {
        // active functions at the evaluation point
        m_bases[0].basis(patchIndex).active_into( u.col(k), actives);
        // for each component, make global DoFs out of these local DoFs
        std::vector< gsMatrix<unsigned> > ci_actives(m_dim,actives);
        for( index_t ci = 0; ci < m_dim; ++ci )
            m_dofMappers[ci].localToGlobal(actives, patchIndex, ci_actives[ci]);

        // compute the gradients of the basis functions at the eval.pt.
        m_bases[0].basis(patchIndex).deriv_into( u.col(k), basisGrad);
        // and transform them to the physical space
        physGrad.resize(m_dim, actives.rows() );
        physGrad.setZero();
        geoEval->evaluateAt( u.col(k) );
        geoEval->transformGradients( index_t(0) , basisGrad, physGrad);

        // compute the partial derivatives of the
        // displacement field
        uPartDers.setZero();
        for( index_t dc=0; dc < m_dim; dc++) // component
            for( index_t dd=0; dd < m_dim; dd++) // derivative
                for( size_t ai=0; ai < size_t( actives.rows() ); ai++ ) // active function
                {
                    unsigned tmpi = ci_actives[ dc ](ai,0);
                    real_t tmpCoef;
                    if( m_dofMappers[dc].is_free_index( tmpi) )
                        tmpCoef = solVector( tmpi, 0);
                    else // dirichlet boundary
                        tmpCoef = m_ddof( m_dofMappers[dc].global_to_bindex(tmpi) );

                    uPartDers( dc, dd ) += tmpCoef * physGrad( dd, ai );
                }

        if( m_dim == 2 )
        {
            // sigma_{11}:
            result(0,k) = m_lambda * ( uPartDers(0,0)+uPartDers(1,1) ) + 2*m_mu*uPartDers(0,0);
            // sigma_{22}:
            result(1,k) = m_lambda * ( uPartDers(0,0)+uPartDers(1,1) ) + 2*m_mu*uPartDers(1,1);
            // sigma_{12}:
            result(2,k) = m_mu * ( uPartDers(0,1) + uPartDers(1,0) );

            if( computeVonMises )
            {
                result(3,k) = math::sqrt(
                            result(0,k)*result(0,k)
                            - result(0,k)*result(1,k)
                            + result(1,k)*result(1,k)
                            + 3.0 * result(2,k)*result(2,k) );
            }

        }
        else if( m_dim == 3 )
        {
            real_t tmp = m_lambda * ( uPartDers(0,0)+uPartDers(1,1)+uPartDers(2,2) );
            result(0,k) = tmp + 2*m_mu * uPartDers(0,0); // sigma_{11}
            result(1,k) = tmp + 2*m_mu * uPartDers(1,1); // sigma_{22}
            result(2,k) = tmp + 2*m_mu * uPartDers(2,2); // sigma_{33}
            result(3,k) = m_mu * ( uPartDers(0,1) + uPartDers(1,0) ); // sigma_{12}
            result(4,k) = m_mu * ( uPartDers(0,2) + uPartDers(2,0) ); // sigma_{13}
            result(5,k) = m_mu * ( uPartDers(2,1) + uPartDers(1,2) ); // sigma_{23}

            if( computeVonMises )
            {
                T tmp2 =   (result(0,k)-result(1,k))*(result(0,k)-result(1,k))
                        + (result(0,k)-result(2,k))*(result(0,k)-result(2,k))
                        + (result(1,k)-result(2,k))*(result(1,k)-result(2,k))
                        + 6.0 * ( result(3,k)*result(3,k) + result(4,k)*result(4,k) + result(5,k)*result(5,k) );

                result(6,k) = math::sqrt( 0.5 * tmp2 );
            }

        }
    }
}

template <class T>
void gsElasticityAssembler<T>::constructStresses(const gsMatrix<T>& solVector,
                                                 gsMultiFunction<T>& result,
                                                 stress_type::type type) const
{
    result.clear();
    for (std::size_t i = 0; i < m_patches.nPatches(); ++i )
    {
        result.addFunction(typename gsFunction<T>::uPtr(new gsStressFunction<T>(solVector,*this,i,type)));
        //gsFunction<T> * temp = new gsStressFunction<T>(solVector,*this,i,type);
        //result.addFunction(temp);
    }
}

template <class T>
const gsMatrix<T> & gsElasticityAssembler<T>::rhs()
{
    if (m_rhsExtra.rows() == 0)
    {
        return m_rhs;
    }
    else
    {
        return m_rhsExtra;
    }
}

template <class T>
void gsElasticityAssembler<T>::addNeummannData(const gsMultiPatch<> &sourceGeometry,
                                               const gsMultiPatch<> &sourceSolution,
                                               int sourcePatch, const boxSide &sourceSide,
                                               int targetPatch, const boxSide &targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceGeometry.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_patches.patch(targetPatch).boundary(targetSide)));
    if (matchingMode != 0)
    {
        if (m_rhsExtra.rows() == 0 && m_rhs.rows() != 0)
            m_rhsExtra = m_rhs;

        const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);
        const gsMatrix<T> coefs = sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs();

        std::map<unsigned,T> neumannDoFs;
        for (int i = 0; i < localBoundaryIndices.rows(); ++i)
        {
            if (matchingMode == 1)
                neumannDoFs[localBoundaryIndices(i)] = coefs(i);
            if (matchingMode == -1)
                neumannDoFs[localBoundaryIndices(i)] = coefs(localBoundaryIndices.rows() - i - 1);
        }

        gsVisitorElasticityPressure<T> visitor(neumannDoFs,targetSide,m_rhsExtra);
        this->apply(visitor,targetPatch,targetSide);
    }
}

template <class T>
void gsElasticityAssembler<T>::addNeummannData(const gsField<> &sourceField,
                                               int sourcePatch, const boxSide &sourceSide,
                                               int targetPatch, const boxSide &targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceField.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_patches.patch(targetPatch).boundary(targetSide)));
    if (matchingMode != 0)
    {
        if (m_rhsExtra.rows() == 0 && m_rhs.rows() != 0)
            m_rhsExtra = m_rhs;

        const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);
        const gsMatrix<T> coefs = sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs();

        std::map<unsigned,T> neumannDoFs;
        for (int i = 0; i < localBoundaryIndices.rows(); ++i)
        {
            if (matchingMode == 1)
                neumannDoFs[localBoundaryIndices(i)] = coefs(i);
            if (matchingMode == -1)
                neumannDoFs[localBoundaryIndices(i)] = coefs(localBoundaryIndices.rows() - i - 1);
        }

        gsVisitorElasticityPressure<T> visitor(neumannDoFs,targetSide,m_rhsExtra);
        this->apply(visitor,targetPatch,targetSide);
    }
}

template <class T>
int gsElasticityAssembler<T>::checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                                      const gsGeometry<> & targetBoundary)
{
    const T absTol = 1e-6; // another magic number

    gsMatrix<> startPoint;
    startPoint.resize(1,1);
    startPoint << 0.;
    gsMatrix<> endPoint;
    endPoint.resize(1,1);
    endPoint << 1.;

    gsMatrix<> sourceStart = sourceBoundary.eval(startPoint);
    gsMatrix<> targetStart = targetBoundary.eval(startPoint);

    gsMatrix<> sourceEnd = sourceBoundary.eval(endPoint);
    gsMatrix<> targetEnd = targetBoundary.eval(endPoint);

    if ((sourceStart-targetStart).norm() + (sourceEnd-targetEnd).norm() < absTol)
    {
        return 1;
    }
    else if ((sourceStart-targetEnd).norm() + (sourceEnd-targetStart).norm() < absTol)
    {
        return -1;
    }
    else
    {
        gsInfo << "Doesn't look like matching boundaries!\n";
        gsInfo << "Source start\n" << sourceStart << std::endl << "Source end\n" << sourceEnd << std::endl;
        gsInfo << "Target start\n" << targetStart << std::endl << "Target end\n" << targetEnd << std::endl;
        return 0;
    }
}

template <class T>
void gsElasticityAssembler<T>::setDirichletDoFs(const gsMultiPatch<> & sourceGeometry,
                                                const gsMultiPatch<> & sourceSolution,
                                                int sourcePatch, const boxSide & sourceSide,
                                                int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceGeometry.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_patches.patch(targetPatch).boundary(targetSide)));
    if (matchingMode == 1)
    {
        setDirichletDoFs(sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs(),targetPatch,targetSide);
    }
    else if (matchingMode == -1)
    {
        setDirichletDoFs(sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs().colwise().reverse(),targetPatch,targetSide);
    }

}

template <class T>
void gsElasticityAssembler<T>::setDirichletDoFs(const gsField<> & sourceField,
                                                int sourcePatch, const boxSide & sourceSide,
                                                int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceField.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_patches.patch(targetPatch).boundary(targetSide)));
    if (matchingMode == 1)
    {
        setDirichletDoFs(sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs(),targetPatch,targetSide);
    }
    else if (matchingMode == -1)
    {
        setDirichletDoFs(sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs().colwise().reverse(),targetPatch,targetSide);
    }

}

template <class T>
void gsElasticityAssembler<T>::setDirichletDoFs(const gsMatrix<> & ddofs,
                                                int targetPatch,
                                                const boxSide & targetSide)
{
    const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);

    for (index_t d = 0; d < m_dim; ++d)
    {
        gsMatrix<unsigned> systemIndices;
        const gsDofMapper & mapper = m_dofMappers[d];
        mapper.localToGlobal(localBoundaryIndices,targetPatch,systemIndices);

        for (index_t i = 0; i != systemIndices.size(); ++i)
        {
            index_t dirichletIndex = mapper.global_to_bindex(systemIndices(i));
            m_ddof(dirichletIndex,0) = ddofs(i,d);
        }
    }
}


template <class T>
void gsElasticityAssembler<T>::clampPatchCorner(int patch,int corner)
{
    unsigned cornerIndex;
    if (corner == 0 || corner == 1)
    {
        const gsMatrix<unsigned> & boundaryIndices = m_bases[0].basis(patch).boundary(boundary::south);
        if (corner == 0)
            cornerIndex = boundaryIndices(0);
        else
            cornerIndex = boundaryIndices(boundaryIndices.rows()-1);
    }
    else if (corner == 2 || corner == 3)
    {
        const gsMatrix<unsigned> & boundaryIndices = m_bases[0].basis(patch).boundary(boundary::north);
        if (corner == 2)
            cornerIndex = boundaryIndices(0);
        else
            cornerIndex = boundaryIndices(boundaryIndices.rows()-1);
    }
    else
    {
        gsInfo << "Bad corner number!\n"
               << "0 - south-west\n"
               << "1 - south-east\n"
               << "2 - north-west\n"
               << "3 - north-east\n";
        return;
    }

    const T PP = 1e9; // magic number

    for (int i = 0; i < m_dim; ++i)
    {
        index_t globalIndex = m_dofMappers[i].index(cornerIndex,patch);
        m_matrix.coeffRef(globalIndex,globalIndex) = PP;
        m_rhs(globalIndex) = 0;
    }
}

template <class T>
void gsElasticityAssembler<T>::deformGeometry(const gsMatrix<T> & solVector,
                                              gsMultiPatch<T> &result)
{
    result.clear();
    gsMatrix<T> newCoeffs;

    for (size_t p = 0; p < m_patches.nPatches(); ++p)
    {
        index_t pNum = m_bases[0][p].size();
        newCoeffs.resize(pNum,m_dim);

        gsMatrix<T> & coeffs = m_patches.patch(p).coefs();

        for (index_t d = 0; d < m_dim; ++d)
        {
            const gsDofMapper & mapper = m_dofMappers[d];

            for (index_t i = 0; i < pNum; ++i )
            {
                if (mapper.is_free(i,p))
                {
                    newCoeffs(i,d) = coeffs(i,d) + solVector(mapper.index(i,p),0);
                }
                else
                {
                    newCoeffs(i,d) = coeffs(i,d) + m_ddof(mapper.bindex(i,p),0);
                }
            }
        }
        result.addPatch(m_bases[0][p].makeGeometry(give(newCoeffs)));
    }
}


template <class T>
void gsStressFunction<T>::eval_into(const gsMatrix< T > & u,gsMatrix< T > & result ) const
{
    gsMatrix<T> fullStress;
    result.clear();
    //gsInfo << m_displacement<<std::endl;
    m_assembler.computeStresses(m_displacement,u,m_patch,fullStress,true);
    result.resize(targetDim(),u.cols());

    switch(m_type)
    {
    case stress_type::normal:
    {
        result.block(0,0,targetDim(),u.cols()) = fullStress.block(0,0,targetDim(),u.cols());
        break;
    }
    case stress_type::shear:
    {
        result.block(0,0,targetDim(),u.cols()) = fullStress.block(domainDim(),0,targetDim(),u.cols());
        break;
    }
    case stress_type::von_mises:
    {
        result.block(0,0,1,u.cols()) = fullStress.block(fullStress.rows()-1,0,1,u.cols());
        break;
    }
    default:
        gsInfo << "Bad stress type\n";
        result.setZero();
        break;
    };
}

}// namespace gismo
