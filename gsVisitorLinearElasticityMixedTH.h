/** @file gsVisitorLinearElasticityMixedTH.h

    @brief Taylor-Hood element visitor for a 2-field mixed method for (near) incompressible linear elasticity for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once


namespace gismo
{


template <class T>
class gsVisitorLinearElasticityMixedTH
{
public:

    /// Constructor
    gsVisitorLinearElasticityMixedTH(T lambda, T mu, T rho, const gsFunction<T> & body_force, T tfac = 1.0) : 
    m_lambda(lambda),
    m_mu(mu),
	m_rho(rho),
    m_bodyForce_ptr(&body_force),
	m_tfac(tfac)
    { }

    void initialize(gsBasisRefs<T> const   & basisRefs,
                    gsQuadRule<T>          & rule, 
                    unsigned               & evFlags )
    {
        m_dim = basisRefs.front().dim();
		m_dimStrain = (m_dim*(m_dim+1))/2;

		gsVector<index_t> numQuadNodes(m_dim);
        for (size_t i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basisRefs.front().degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_GRAD_TRANSFORM;

		m_C.resize(m_dimStrain,m_dimStrain);
		if (m_dim == 2)
        {
			m_C(0,0) = m_C(1,1) = 2*m_mu;
			m_C(2,2) = m_mu;
		}
		else if (m_dim == 3)
		{
			m_C(0,0) = m_C(1,1) = m_C(2,2) = 2*m_mu;
			m_C(3,3) = m_C(4,4) = m_C(5,5) = m_mu;
		}
    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basisRefs.front().active_into(quNodes.col(0), actives);
        numActive   = actives.rows();
		basisRefs.back(). active_into(quNodes.col(0), actives_p);
        numActive_p = actives_p.rows();
        
        // Evaluate basis functions on element
        basisRefs.front().evalAllDers_into( quNodes, 1, basisData);
		basisRefs.back(). eval_into( quNodes, basisVals_p);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);
        
        // Evaluate right-hand side at the geometry points
        m_bodyForce_ptr->eval_into( geoEval.values(), forceVals );
        
        // Initialize local matrix/rhs
        localMatK.setZero(m_dim*numActive, m_dim*numActive);
		localMatB.setZero(numActive_p, m_dim*numActive);
		localMatC.setZero(numActive_p, numActive_p);
        localRhs_u.setZero(m_dim*numActive, 1);		
        localRhs_p.setZero(numActive_p, 1);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const typename gsMatrix<T>::Block bVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block bGrads = basisData.middleRows(numActive, m_dim*numActive);

		const typename gsMatrix<T>::Block bVals_p  = basisVals_p.topRows(numActive_p);

		const T v_2mulam = 2.*m_mu;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {           
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);

			// compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);

			if (m_dim == 2)
			{
				for (index_t i = 0; i < numActive; i++)
				{
					for (index_t j = 0; j < numActive; j++)
					{
						localMatK(0*numActive+i, 0*numActive+j) += weight * 
							( v_2mulam*physGrad(0,i)*physGrad(0,j) + m_mu*physGrad(1,i)*physGrad(1,j) );
						localMatK(0*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(0,j) );
						localMatK(1*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(1,j) );
						localMatK(1*numActive+i, 1*numActive+j) += weight * 
							( v_2mulam*physGrad(1,i)*physGrad(1,j) + m_mu*physGrad(0,i)*physGrad(0,j) );
					}

					for (index_t j = 0; j < numActive_p; j++)
					{
						localMatB(j, 0*numActive+i) += weight * m_mu * physGrad(0,i) * basisVals_p(j,k);
						localMatB(j, 1*numActive+i) += weight * m_mu * physGrad(1,i) * basisVals_p(j,k);
					}
				}
			}
			else if (m_dim == 3)
			{
				for (index_t i = 0; i < numActive; i++)
				{
					for (index_t j = 0; j < numActive; j++)
					{
						localMatK(0*numActive+i, 0*numActive+j) += weight * 
							( v_2mulam*physGrad(0,i)*physGrad(0,j) + m_mu*(physGrad(1,i)*physGrad(1,j)+physGrad(2,i)*physGrad(2,j)) );
						localMatK(0*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(0,j) );
						localMatK(0*numActive+i, 2*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(2,j) + m_mu*physGrad(2,i)*physGrad(0,j) );
						localMatK(1*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(1,j) );
						localMatK(1*numActive+i, 1*numActive+j) += weight * 
							( v_2mulam*physGrad(1,i)*physGrad(1,j) + m_mu*(physGrad(0,i)*physGrad(0,j)+physGrad(2,i)*physGrad(2,j)) );
						localMatK(1*numActive+i, 2*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(2,j) + m_mu*physGrad(2,i)*physGrad(1,j) );
						localMatK(2*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(2,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(2,j) );
						localMatK(2*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(2,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(2,j) );
						localMatK(2*numActive+i, 2*numActive+j) += weight * 
							( v_2mulam*physGrad(2,i)*physGrad(2,j) + m_mu*(physGrad(1,i)*physGrad(1,j)+physGrad(0,i)*physGrad(0,j)) );				
					}

					for (index_t j = 0; j < numActive_p; j++)
					{
						localMatB(j, 0*numActive+i) += weight * m_mu * physGrad(0,i) * basisVals_p(j,k);
						localMatB(j, 1*numActive+i) += weight * m_mu * physGrad(1,i) * basisVals_p(j,k);
						localMatB(j, 2*numActive+i) += weight * m_mu * physGrad(2,i) * basisVals_p(j,k);
					}
				}
			}

			// near incompressible
			if (m_lambda < std::numeric_limits<T>::infinity())
			{
				for (index_t i = 0; i < numActive_p; i++)
				{
					for (index_t j = 0; j < numActive_p; j++)
					{
						localMatC(i, j) -= weight * m_mu*m_mu/m_lambda * basisVals_p(i,k) * basisVals_p(j,k);
					}
				}
			}
			
			// Local rhs vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs_u.middleRows(j*numActive,numActive).noalias() += 
                    weight * m_rho * forceVals(j,k) * m_tfac * bVals.col(k) ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }
    
    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const int patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
       	
		// Local DoFs to global DoFs
		std::vector< gsMatrix<unsigned> > ci_actives(m_dim+1,actives);

        for (size_t ci = 0; ci != m_dim; ++ci)
			mappers[ci].localToGlobal(actives, patchIndex, ci_actives[ci]);
		
		mappers[m_dim].localToGlobal(actives_p, patchIndex, ci_actives[m_dim]);

        for (size_t ci = 0; ci!= m_dim; ++ci)
		{          
			for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t gi = ci * numActive +  ai; // row index
                const index_t ii = ci_actives[ci](ai);

                if ( mappers[ci].is_free_index(ii) )
                {
                    rhsMatrix.row(ii) += localRhs_u.row(gi);
                    
					// matrix A
                    for (size_t cj = 0; cj!= m_dim; ++cj)
					{
                        for (index_t aj=0; aj < numActive; ++aj)
                        {
                            const index_t gj = cj * numActive +  aj; // column index
                            const index_t jj = ci_actives[cj](aj);
                            
                            if ( mappers[cj].is_free_index(jj) )
                            {
                                sysMatrix.coeffRef(ii, jj) += localMatK(gi, gj);
                            }
                            else // Fixed DoF ?
                            {
                                const index_t bjj = mappers[cj].global_to_bindex(jj);
								rhsMatrix.row(ii).noalias() -= localMatK(gi, gj) * 
                                    eliminatedDofs.row( bjj );
                            }
                        }
					}

					// matrix B and B'
					size_t cj = m_dim;
                    for (index_t aj=0; aj < numActive_p; ++aj)
                    {
                        const index_t gj = aj;					// column index
                        const index_t jj = ci_actives[cj](aj);
                            
                        if ( mappers[cj].is_free_index(jj) )
                        {
                            sysMatrix.coeffRef(ii, jj) += localMatB(gj, gi);
							sysMatrix.coeffRef(jj, ii) += localMatB(gj, gi);
                        }
                        else // Fixed DoF ?
                        {
                            const index_t bjj = mappers[cj].global_to_bindex(jj);
							rhsMatrix.row(ii).noalias() -= localMatB(gj, gi) * 
                                eliminatedDofs.row( bjj );
                        }
                    }
					
                }
            }

		}

		// matrix C
		size_t ci = m_dim;
		for (index_t ai=0; ai < numActive_p; ++ai)
        {
            const index_t gi = ai;					// row index
            const index_t ii = ci_actives[ci](ai);

			if ( mappers[ci].is_free_index(ii) )
            {
                rhsMatrix.row(ii) += localRhs_p.row(gi);

				size_t cj = m_dim;
                for (index_t aj=0; aj < numActive_p; ++aj)
                {
                    const index_t gj = aj;					// column index
                    const index_t jj = ci_actives[cj](aj);
                            
                    if ( mappers[cj].is_free_index(jj) )
                    {
                        sysMatrix.coeffRef(ii, jj) += localMatC(gi, gj);
                    }                    
                }
			}
		}

    }
    

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    // Basis values
    gsMatrix<T>        basisData;
    gsMatrix<unsigned> actives;
	gsMatrix<T>		   physGrad, physGrad_symm;
    index_t            numActive;

	// Pressure basis values
    gsMatrix<T>        basisVals_p;
    gsMatrix<unsigned> actives_p;
	index_t            numActive_p;

    gsVector<T> normal;

    // Material matrix
    gsMatrix<T> m_C;
    

protected:

	// Dimension
	size_t m_dim, m_dimStrain;

    // Lambda, mu, rho
    T m_lambda, m_mu, m_rho;

	// Factor for time-dependent body force
	T m_tfac;

protected:

    // Surface forces
    const gsFunction<T> * m_bodyForce_ptr;

    // Local values of the surface forces
    gsMatrix<T> forceVals;
    
protected:
    // Local matrices
    gsMatrix<T> localMatK, localMatB, localMatC;
    gsMatrix<T> localRhs_u, localRhs_p;
};


} // namespace gismo

