/** @file gsVisitorNonLinElasticity.h

    @brief Element visitor for nonlinear elasticity for 2D plain strain 
	and 3D continua.

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
                              const gsGeometry<T> & deformed,
							  T tfac = 1.0) : 
    Base(lambda,mu,rho,body_force,tfac),
    m_deformation(deformed.evaluator(NEED_JACOBIAN)) // NEED_MEASURE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM))
    { 
		m_dim = body_force.targetDim();

		m_MATERIAL_LAW = 1;

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
        m_deformation = memory::make_unique(deformed.evaluator(NEED_JACOBIAN)); // (NEED_MEASURE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM));
    }

	/// Sets the \em m_MATERIAL_LAW to 0: St. Venant-Kirchhoff, 1: Neo-Hooke
    void set_MaterialLaw(const int material)
    {
        GISMO_ASSERT( (material >= 0) && (material <= 1), "Only 0 (St. Venant-Kirchhoff) or 1 (Neo-Hooke) allowed!");
		m_MATERIAL_LAW = material;
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
		if (m_MATERIAL_LAW == 0)
			assemble_StVK(element, geoEval, quWeights);
		else
			assemble_NeoH(element, geoEval, quWeights);
	}

    inline void assemble_NeoH(gsDomainIterator<T>    & element, 
                             gsGeometryEvaluator<T> & geoEval,
                             gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];
		//const gsMatrix<T> & defGrads = m_deformation->jacobians();

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
                        logdetF = math::log(detF);
			mulamlogJ = m_mu - m_lambda * logdetF;

			// Inverse of Fi = F^-1
			defGrad_inv = defGrad.inverse();

			// Local internal force vector contribution, mu*(F-Fi') + lambda*log(J)*Fi' 
			locResMat = (weight * m_mu) * (defGrad - defGrad_inv.transpose()) + (weight * m_lambda * logdetF) * defGrad_inv.transpose();

			// 1st basis function (U/i)
			for (index_t i = 0; i < numActive; i++)
			{		
				locResVec = locResMat * physGrad.col(i);

				// Spatial dimensions of 1st basis function
				for (size_t di = 0; di < m_dim; di++)
				{				
					// Write to Rhs
					localRhs(di*numActive+i) -= locResVec(di);
					
					// Write gradient as matrix
					gradU.setZero();
					gradU.row(di) = physGrad.col(i);

					defGrad_inv_gradU = defGrad_inv * gradU;				// Fi*dU
					defGrad_inv_gradU_trace = defGrad_inv_gradU.trace();	// tr(Fi*dU) = Fi':dU
						
					// 2nd basis function (V/j)
					//for (index_t j = 0; j < numActive; j++)
					// Exploit symmetry of K
					for (index_t j = i; j < numActive; j++)
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
                    weight * m_rho * forceVals(j,k) * m_tfac * bVals.col(k) ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }

	inline void assemble_StVK(gsDomainIterator<T>    & element, 
                             gsGeometryEvaluator<T> & geoEval,
                             gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];
		//const gsMatrix<T> & defGrads = m_deformation->jacobians();

		size_t v_dim = m_dim*(m_dim+1)/2;

		gsMatrix<T>	strainRCG(m_dim,m_dim), strainGLG(m_dim,m_dim), stressPK2(m_dim,m_dim);
		gsVector<T>	vec_RCG(v_dim), vec_GLG(v_dim), vec_PK2(v_dim);
		gsMatrix<T> mat_D(v_dim,v_dim);
		gsMatrix<T> mat_B1(v_dim,m_dim), mat_B2(v_dim,m_dim);
		gsVector<T> vec_geo(m_dim);
		gsMatrix<T> Ktg_mat1(m_dim,v_dim);
		T fac_geo;
		gsMatrix<T>	locKtgMat(m_dim,m_dim);

		mat_D.setZero();

		if (m_dim == 2)
		{	
			mat_D(0,0) = mat_D(1,1) = 2*m_mu + m_lambda;
			mat_D(0,1) = mat_D(1,0) = m_lambda;
			mat_D(2,2) = m_mu;
		}
		else if (m_dim == 3)
		{			
			mat_D(0,0) = mat_D(1,1) = mat_D(2,2) = 2*m_mu + m_lambda;
			mat_D(0,1) = mat_D(0,2) = mat_D(1,2) = m_lambda;
			mat_D(1,0) = mat_D(2,0) = mat_D(2,1) = m_lambda;
			mat_D(3,3) = mat_D(4,4) = mat_D(5,5) = m_mu;
		}

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

			// Right Cauchy Green strain, C = F'*F
			strainRCG = defGrad.transpose() * defGrad;

			// Green-Lagrange Green strain, E = 0.5*(C-I)
			strainGLG = strainRCG;
			for (size_t di = 0; di < m_dim; di++)
				strainGLG(di,di) -= 1.;
			strainGLG *= 0.5;

			// Green-Lagrange Green strain as vector
			if (m_dim == 2)
			{
				vec_GLG(0) = strainGLG(0,0);
				vec_GLG(1) = strainGLG(1,1);
				vec_GLG(2) = 2.*strainGLG(0,1);
				
				vec_PK2 = mat_D * vec_GLG;

				stressPK2(0,0) = vec_PK2(0);
				stressPK2(1,1) = vec_PK2(1);
				stressPK2(0,1) = vec_PK2(2);
				stressPK2(1,0) = vec_PK2(2);
			}
			else if (m_dim == 3)
			{
				vec_GLG(0) = strainGLG(0,0);
				vec_GLG(1) = strainGLG(1,1);
				vec_GLG(2) = strainGLG(2,2);
				vec_GLG(3) = 2.*strainGLG(0,1);
				vec_GLG(4) = 2.*strainGLG(1,2);
				vec_GLG(5) = 2.*strainGLG(0,2);

				vec_PK2 = mat_D * vec_GLG;

				stressPK2(0,0) = vec_PK2(0);
				stressPK2(1,1) = vec_PK2(1);
				stressPK2(2,2) = vec_PK2(2);
				stressPK2(0,1) = vec_PK2(3);
				stressPK2(1,2) = vec_PK2(4);
				stressPK2(0,2) = vec_PK2(5);
				stressPK2(1,0) = vec_PK2(3);
				stressPK2(2,1) = vec_PK2(4);
				stressPK2(2,0) = vec_PK2(5);
			}

			// 1st basis function (U/i)
			for (index_t i = 0; i < numActive; i++)
			{		
				Set_mat_B(mat_B1, defGrad, physGrad.col(i));

				locResVec = mat_B1.transpose() * vec_PK2;

				for (size_t di = 0; di < m_dim; di++)
					localRhs(di*numActive+i) -= weight * locResVec(di);

				// Geometric tangent 					
				vec_geo = stressPK2 * physGrad.col(i);

				// Material tangent				
				Ktg_mat1 = mat_B1.transpose() * mat_D;

				// 2nd basis function (V/j)
				//for (index_t j = 0; j < numActive; j++)
				// Exploit symmetry of K
				for (index_t j = i; j < numActive; j++)
				{						
					Set_mat_B(mat_B2, defGrad, physGrad.col(j));
					
					// Material tangent		
					locKtgMat = Ktg_mat1 * mat_B2;

					// Geometric tangent 					
					fac_geo = vec_geo.transpose() * physGrad.col(j);

					for (size_t dj = 0; dj < m_dim; dj++)
						locKtgMat(dj,dj) += fac_geo;

					// Write to local matrix 
					for (size_t di = 0; di < m_dim; di++)
						for (size_t dj = 0; dj < m_dim; dj++)
							localMat(di*numActive+i, dj*numActive+j) += weight * locKtgMat(di,dj);
				}
			
			}

			// Local external force vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                    weight * m_rho * forceVals(j,k) * m_tfac * bVals.col(k) ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }


	// We need to over-write localToGlobal from linear visitor, 
	//   because RHS should not be modified anymore for fixed DoFs!
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

        for (index_t ai=0; ai < numActive; ++ai)
		{          
			for (size_t ci = 0; ci!= m_dim; ++ci)
            {
                const index_t gi = ci * numActive +  ai; // row index
                const index_t ii = ci_actives[ci](ai);

                if ( mappers[ci].is_free_index(ii) )
                {
                    rhsMatrix.row(ii) += localRhs.row(gi);
                    
					//for (index_t aj=0; aj < numActive; ++aj)
					// Exploit symmetry of K
					for (index_t aj=ai; aj < numActive; ++aj)
					{
						for (size_t cj = 0; cj!= m_dim; ++cj)                        
                        {
                            const index_t gj = cj * numActive +  aj; // column index
                            const index_t jj = ci_actives[cj](aj);
                            
                            if ( mappers[cj].is_free_index(jj) )
                            {
                                sysMatrix.coeffRef(ii, jj) += localMat(gi, gj);
								if (aj > ai)
									sysMatrix.coeffRef(jj, ii) += localMat(gi, gj);
                            }
                        }
					}
                }
				
            }

		}
    }

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:	
	inline void Set_mat_B( gsMatrix<T> &B, const gsMatrix<T> &dg, const gsVector<T> &mg )
	{
		if (m_dim == 2)
		{
			B(0,0) = dg(0,0) * mg(0);	
			B(0,1) = dg(1,0) * mg(0);	
			B(1,0) = dg(0,1) * mg(1);	
			B(1,1) = dg(1,1) * mg(1);	
			B(2,0) = dg(0,0) * mg(1) + dg(0,1) * mg(0);
			B(2,1) = dg(1,0) * mg(1) + dg(1,1) * mg(0);
		}
		else if (m_dim == 3) 
		{
			for (size_t i = 0; i < m_dim; i++)
			{
				size_t k = (i+1)%m_dim;
				for (size_t j = 0; j < m_dim; j++)
				{
					B(i,j) = dg(j,i) * mg(i);					
					B(i+m_dim,j) = dg(j,i) * mg(k) + dg(j,k) * mg(i);
				}
			}
		}
	}

protected:
	// Lambda, mu, rho
    using Base::m_lambda;
	using Base::m_mu;
	using Base::m_rho;

	// Body forces
    using Base::m_bodyForce_ptr;

	/// Contains the geometry evaluations for the deformed configuration
    typename gsGeometry<T>::Evaluator m_deformation;

	// Body forces
	using Base::m_tfac;

	// Type of material law (0: St. Venant-Kirchhoff, 1: Neo-Hooke)
	int m_MATERIAL_LAW;

protected:
	// Dimension
	using Base::m_dim;
	using Base::m_dimStrain;
	
	// Basis values
    using Base::basisData;
    using Base::actives;
    using Base::numActive;
	using Base::physGrad;
	
    // Local values of the surface forces
    using Base::forceVals;

protected:
    // Local matrices
    using Base::localMat;
    using Base::localRhs;	

protected:
	// Kinematics
	T weight;
	gsMatrix<T> defDer_k;
	gsVector<T> defDer_kv;

	gsMatrix<T> displGrad;
	gsMatrix<T> defGrad;
	gsMatrix<T> defGrad_inv;
	T detF, logdetF, mulamlogJ;

	gsMatrix<T> gradU;
	gsMatrix<T> gradV;
	gsMatrix<T> defGrad_inv_gradU;
	gsMatrix<T> defGrad_inv_gradV;
	T defGrad_inv_gradU_trace;

	gsMatrix<T> locResMat;
	gsVector<T> locResVec;
	T locKtgVal;
};


} // namespace gismo

