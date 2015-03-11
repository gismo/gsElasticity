/** @file gsVisitorNonLinElasticity.h

    @brief Element visitor for nonlinear elasticity for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once

#include <gsElasticity/gsVisitorLinearElasticity.h>

namespace gismo
{


template <class T>
class gsVisitorNonLinElasticity : public gsVisitorLinearElasticity<T>
{
public:

	typedef gsVisitorLinearElasticity<T> Base;

    /// Constructor
    gsVisitorNonLinElasticity(T lambda, T mu, T rho, 
		                      const gsFunction<T> & body_force, 
                              const gsGeometry<T> & deformed) : 
    Base(lambda,mu,rho,body_force),
    m_deformation(deformed.evaluator(NEED_MEASURE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM))
    { 
		m_dim = body_force.targetDim();

		defDer_k.resize(m_dim,m_dim);
		defDer_kv.resize(m_dim*m_dim);

		displGrad.resize(m_dim,m_dim);
		defGrad.resize(m_dim,m_dim);
		defGrad_inv.resize(m_dim,m_dim);

		gradU.resize(m_dim,m_dim);
		gradV.resize(m_dim,m_dim);
		defGrad_inv_gradU.resize(m_dim,m_dim);
		defGrad_inv_gradV.resize(m_dim,m_dim);
		
		locResMat.resize(m_dim,m_dim);
		locResVec.resize(m_dim);
	}

	/// Sets the gsGeometryEvaluator \em m_deformation using \em deformed
    void setDeformed(const gsGeometry<T> & deformed)
    {
        m_deformation = safe(deformed.evaluator(NEED_MEASURE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM));
    }

   /// Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        Base::evaluate(basis,geoEval,quNodes);

        // Evaluate deformed shell 
        m_deformation->evaluateAt(quNodes);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const typename gsMatrix<T>::Block bVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block bGrads = basisData.middleRows(numActive, m_dim*numActive);
		const gsMatrix<T> defGrads = m_deformation->jacobians();

		/*
		T weight;
		gsMatrix<T> defDer_k(m_dim,m_dim);
		gsVector<T> defDer_kv(m_dim*m_dim);

		gsMatrix<T> displGrad(m_dim,m_dim);
		gsMatrix<T> defGrad(m_dim,m_dim);
		gsMatrix<T> defGrad_inv(m_dim,m_dim);
		T detF, mulamlogJ;

		gsMatrix<T> gradU(m_dim,m_dim);
		gsMatrix<T> gradV(m_dim,m_dim);
		gsMatrix<T> defGrad_inv_gradU(m_dim,m_dim);
		gsMatrix<T> defGrad_inv_gradV(m_dim,m_dim);
		T defGrad_inv_gradU_trace;

		gsMatrix<T> locResMat(m_dim,m_dim);
		gsVector<T> locResVec(m_dim);
		T locKtgVal;
		*/

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {           
            // Multiply weight by the geometry measure
            weight = quWeights[k] * geoEval.measure(k);

			// compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);

			//geoEval.transformGradients(k, 
			//gsAsConstMatrix<T,3,3> geoJac(geoEval->jacobian(k).data(),3,3 );
			//gsAsConstMatrix<T,3,3> defDer(m_deformation->jacobian(k).data(),3,3 );
			defDer_k = m_deformation->jacobian(k).transpose();
			//defDer_kv = defDer_k.asVector();
			//geoEval.transformGradients(k, defDer_kv, displGrad);

			// Displacement gradient, H = du = (dx/dxi)'^-1 * (du/dxi)'
			displGrad = geoEval.jacobian(k).transpose().inverse() * defDer_k;

			// Deformation gradient, F = I + H
			defGrad = displGrad.transpose();
			for (size_t di = 0; di < m_dim; di++)
				defGrad(di,di) += 1.;

			// Determinant of deformation gradient, J = det(F)
			detF = defGrad.determinant();
			mulamlogJ = m_mu - m_lambda * std::log(detF);

			// Inverse of Fi = F^-1
			defGrad_inv = defGrad.inverse();

			// 1st basis function (U/i)
			for (index_t i = 0; i < numActive; i++)
			{
				// Local internal force vector contribution, mu*(F-Fi') + lambda*log(J)*Fi' 
				locResMat = m_mu * (defGrad - defGrad_inv.transpose()) + m_lambda * std::log(detF) * defGrad_inv.transpose();
				locResVec = locResMat * physGrad.col(i);

				// Spatial dimensions of 1st basis function
				for (size_t di = 0; di < m_dim; di++)
				{				
					// Write to Rhs
					localRhs(di*numActive+i) -= weight * locResVec(di);
					
					// Write gradient as matrix
					gradU.setZero();
					gradU.row(di) = physGrad.col(i);

					defGrad_inv_gradU = defGrad_inv * gradU;				// Fi*dU
					defGrad_inv_gradU_trace = defGrad_inv_gradU.trace();	// tr(Fi*dU) = Fi':dU
						
					// 2nd basis function (V/j)
					for (index_t j = 0; j < numActive; j++)
					{
						
						// Spatial dimensions of 2nd basis function
						for (size_t dj = 0; dj < m_dim; dj++)
						{				
							// Write gradient as matrix
							gradV.setZero();
							gradV.row(dj) = physGrad.col(j);

							defGrad_inv_gradV = defGrad_inv * gradV;		// Fi*dV

							// Local tangent stiffnees matrix contribution
							locKtgVal = m_mu * ( gradU.transpose()*gradV ).trace()
							          + mulamlogJ * ( defGrad_inv_gradU * defGrad_inv_gradV ).trace()
							          + m_lambda * defGrad_inv_gradU_trace * defGrad_inv_gradV.trace();

							// Write to Mat
							localMat(di*numActive+i, dj*numActive+j) += weight * locKtgVal;
						}
					}
				}				
			}

			// Local external force vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                    weight * m_rho * forceVals(j,k) * bVals.col(k) ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }
    
	/*
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
    */

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

	/// Contains the geometry evaluations for the deformed configuration
    typename gsGeometry<T>::Evaluator m_deformation;
	
	// Kinematics
	T weight;
	gsMatrix<T> defDer_k;
	gsVector<T> defDer_kv;

	gsMatrix<T> displGrad;
	gsMatrix<T> defGrad;
	gsMatrix<T> defGrad_inv;
	T detF, mulamlogJ;

	gsMatrix<T> gradU;
	gsMatrix<T> gradV;
	gsMatrix<T> defGrad_inv_gradU;
	gsMatrix<T> defGrad_inv_gradV;
	T defGrad_inv_gradU_trace;

	gsMatrix<T> locResMat;
	gsVector<T> locResVec;
	T locKtgVal;

protected:

    // Basis values
    using Base::basisData;
    using Base::actives;
    using Base::numActive;
    using Base::normal;
	using Base::physGrad;	

protected:

	// Dimension
	using Base::m_dim;
	using Base::m_dimStrain;

    // Lambda, mu, rho
    using Base::m_lambda;
	using Base::m_mu;
	using Base::m_rho;

protected:

    // Surface forces
    using Base::m_bodyForce_ptr;

    // Local values of the surface forces
    using Base::forceVals;
    
protected:
    // Local matrices
    using Base::localMat;
    using Base::localRhs;

};


} // namespace gismo

