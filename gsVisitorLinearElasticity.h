/** @file gsVisitorLinearElasticity.h

    @brief Element visitor for linear elasticity for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once

//#include <gsElasticity/gsElasticityAssembler.hpp>

namespace gismo
{


template <class T>
class gsVisitorLinearElasticity
{
public:

    /// Constructor
    gsVisitorLinearElasticity(T lambda, T mu, T rho, const gsFunction<T> & body_force, T tfac = 1.0) : 
    m_lambda(lambda),
    m_mu(mu),
	m_rho(rho),
    m_bodyForce_ptr(&body_force)
    {
	   m_tfac = tfac;
	}

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        m_dim = basis.dim();
		m_dimStrain = (m_dim*(m_dim+1))/2;

		gsVector<index_t> numQuadNodes(m_dim);
        for (size_t i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_GRAD_TRANSFORM;

		m_C.resize(m_dimStrain,m_dimStrain);
		if (m_dim == 2)
        {
			m_C(0,0) = m_C(1,1) = 2*m_mu + m_lambda;
			m_C(0,1) = m_C(1,0) = m_lambda;
			m_C(2,2) = m_mu;
		}
		else if (m_dim == 3)
		{
			m_C(0,0) = m_C(1,1) = m_C(2,2) = 2*m_mu + m_lambda;
			m_C(0,1) = m_C(0,2) = m_C(1,2) = m_lambda;
			m_C(1,0) = m_C(2,0) = m_C(2,1) = m_lambda;
			m_C(3,3) = m_C(4,4) = m_C(5,5) = m_mu;
		}
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);
        numActive = actives.rows();
        
        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 1, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);
        
        // Evaluate right-hand side at the geometry points
        m_bodyForce_ptr->eval_into( geoEval.values(), forceVals );
        
        // Initialize local matrix/rhs
        localMat.setZero(m_dim*numActive, m_dim*numActive);
        localRhs.setZero(m_dim*numActive, 1);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const typename gsMatrix<T>::Block bVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block bGrads = basisData.middleRows(numActive, m_dim*numActive);

		const T v_2mulam = 2.*m_mu + m_lambda;

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
						localMat(0*numActive+i, 0*numActive+j) += weight * 
							( v_2mulam*physGrad(0,i)*physGrad(0,j) + m_mu*physGrad(1,i)*physGrad(1,j) );
						localMat(0*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(0,j) );
						localMat(1*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(1,j) );
						localMat(1*numActive+i, 1*numActive+j) += weight * 
							( v_2mulam*physGrad(1,i)*physGrad(1,j) + m_mu*physGrad(0,i)*physGrad(0,j) );
					}
				}
			}
			else if (m_dim == 3)
			{
				for (index_t i = 0; i < numActive; i++)
				{
					for (index_t j = 0; j < numActive; j++)
					{
						localMat(0*numActive+i, 0*numActive+j) += weight * 
							( v_2mulam*physGrad(0,i)*physGrad(0,j) + m_mu*(physGrad(1,i)*physGrad(1,j)+physGrad(2,i)*physGrad(2,j)) );
						localMat(0*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(0,j) );
						localMat(0*numActive+i, 2*numActive+j) += weight * 
							( m_lambda*physGrad(0,i)*physGrad(2,j) + m_mu*physGrad(2,i)*physGrad(0,j) );
						localMat(1*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(1,j) );
						localMat(1*numActive+i, 1*numActive+j) += weight * 
							( v_2mulam*physGrad(1,i)*physGrad(1,j) + m_mu*(physGrad(0,i)*physGrad(0,j)+physGrad(2,i)*physGrad(2,j)) );
						localMat(1*numActive+i, 2*numActive+j) += weight * 
							( m_lambda*physGrad(1,i)*physGrad(2,j) + m_mu*physGrad(2,i)*physGrad(1,j) );
						localMat(2*numActive+i, 0*numActive+j) += weight * 
							( m_lambda*physGrad(2,i)*physGrad(0,j) + m_mu*physGrad(0,i)*physGrad(2,j) );
						localMat(2*numActive+i, 1*numActive+j) += weight * 
							( m_lambda*physGrad(2,i)*physGrad(1,j) + m_mu*physGrad(1,i)*physGrad(2,j) );
						localMat(2*numActive+i, 2*numActive+j) += weight * 
							( v_2mulam*physGrad(2,i)*physGrad(2,j) + m_mu*(physGrad(1,i)*physGrad(1,j)+physGrad(0,i)*physGrad(0,j)) );				
					}
				}
			}

			/*
			if (m_dim == 2)
			{
				m_virtualStrain.row(0) = physGrad.row(0);
				m_virtualStrain.row(1) = physGrad.row(1);
				m_virtualStrain.row(2) = physGrad.row(0) + physGrad.row(1);
			}
			else if (m_dim == 3)
			{
				m_virtualStrain.row(0) = physGrad.row(0);
				m_virtualStrain.row(1) = physGrad.row(1);
				m_virtualStrain.row(2) = physGrad.row(2);
				m_virtualStrain.row(3) = physGrad.row(0) + physGrad.row(1);
				m_virtualStrain.row(4) = physGrad.row(1) + physGrad.row(2);
				m_virtualStrain.row(5) = physGrad.row(0) + physGrad.row(2);
			}

			// Local block A
            localMat..noalias()  += weight * (m_virtualStrain.transpose() * m_C * m_virtualStrain);

			// Compute needed quantities...
            computeMaterialMatrix(geoEval, k);
            computeStrainDers(geoEval, bGrads, k);

			// Local stiffness matrix contribution
            localMat.noalias() += weight *  ( E_m_der.transpose() * m_C * E_m_der + 
                                              E_f_der.transpose() * m_C * E_f_der * 
                                              (m_thickness*m_thickness/3.0) );
            */
			
			// Local rhs vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
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
		std::vector< gsMatrix<unsigned> > ci_actives(m_dim,actives);

        for (size_t ci = 0; ci != m_dim; ++ci)
			mappers[ci].localToGlobal(actives, patchIndex, ci_actives[ci]);

        for (size_t ci = 0; ci!= m_dim; ++ci)
		{          
			for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t gi = ci * numActive +  ai; // row index
                const index_t ii = ci_actives[ci](ai);

                if ( mappers[ci].is_free_index(ii) )
                {
                    rhsMatrix.row(ii) += localRhs.row(gi);
                    
                    for (size_t cj = 0; cj!= m_dim; ++cj)
                        for (index_t aj=0; aj < numActive; ++aj)
                        {
                            const index_t gj = cj * numActive +  aj; // column index
                            const index_t jj = ci_actives[cj](aj);
                            
                            if ( mappers[cj].is_free_index(jj) )
                            {
                                sysMatrix.coeffRef(ii, jj) += localMat(gi, gj);
                            }
                            else // Fixed DoF ?
                            {
                                const index_t bjj = mappers[cj].global_to_bindex(jj);
								rhsMatrix.row(ii).noalias() -= localMat(gi, gj) * 
                                    eliminatedDofs.row( bjj );
                            }
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
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
};


} // namespace gismo

