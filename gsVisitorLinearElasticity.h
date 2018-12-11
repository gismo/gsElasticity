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
#include <gsCore/gsFuncData.h>

namespace gismo
{


template <class T>
class gsVisitorLinearElasticity
{
public:

    /// Constructor
    gsVisitorLinearElasticity(T lambda, T mu, T rho, const gsFunction<T> & body_force, T tfac = 1.0,
                              gsFunction<T> * E = nullptr, gsFunction<T> * pr = nullptr) :
    m_lambda(lambda),
    m_mu(mu),
	m_rho(rho),
    m_bodyForce_ptr(&body_force),
    m_tfac(tfac),
    f_E(E),
    f_pr(pr)
    { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        m_dim = basis.dim();
		m_dimStrain = (m_dim*(m_dim+1))/2;

		gsVector<index_t> numQuadNodes(m_dim);
        for (size_t i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometry<T> const & geo,
                         gsMatrix<T> const      & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);
        numActive = actives.rows();

        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 1, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        if (f_E != nullptr)
            f_E->eval_into(quNodes,eData);
        if (f_pr != nullptr)
            f_pr->eval_into(quNodes,prData);

        // Evaluate right-hand side at the geometry points
        m_bodyForce_ptr->eval_into( md.values[0], forceVals );

        // Initialize local matrix/rhs
        localMat.setZero(m_dim*numActive, m_dim*numActive);
        localRhs.setZero(m_dim*numActive, 1);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        T v_2mulam = 2.*m_mu + m_lambda;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            if (f_E != nullptr)
            {
                T E = eData.at(k);
                T pr = prData.at(k);

                m_lambda = E * pr / ( (1. + pr) * (1. - 2.*pr) );
                m_mu = E / (2.*(1.+pr));

                v_2mulam = 2.*m_mu + m_lambda;
            }


            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);

			if (m_dim == 2)
			{
				for (index_t i = 0; i < numActive; i++)
				{
					//for (index_t j = 0; j < numActive; j++)
					// Exploit symmetry of K
					for (index_t j = i; j < numActive; j++)
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
					//for (index_t j = 0; j < numActive; j++)
					// Exploit symmetry of K
					for (index_t j = i; j < numActive; j++)
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

			// Local rhs vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() +=
                    weight * m_rho * forceVals(j,k) * m_tfac * bVals.col(k) ;
        }
       // gsInfo<< "local Mat: \n"<< localMat << "\n";
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        std::vector< gsMatrix<unsigned> > ci_actives(m_dim,actives);
        for (size_t ci = 0; ci != m_dim; ++ci)
            system.mapColIndices(actives, patchIndex, ci_actives[ci], ci);

        system.push(localMat, localRhs, ci_actives,eliminatedDofs);
    }


    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    // Lambda, mu, rho
    T m_lambda, m_mu, m_rho;

    // Body forces
    const gsFunction<T> * m_bodyForce_ptr;

	// Factor for time-dependent body force
	T m_tfac;

protected:
	// Dimension
	size_t m_dim, m_dimStrain;

    // Basis values
    std::vector<gsMatrix<T> > basisData;

    gsMatrix<T> eData;
    gsMatrix<T> prData;

    gsMatrix<unsigned> actives;
	gsMatrix<T>		   physGrad, physGrad_symm;
    index_t            numActive;

    // Material matrix
    gsMatrix<T> m_C;

    // Local values of the surface forces
    gsMatrix<T> forceVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsFunction<T> * f_E;
    gsFunction<T> * f_pr;

    gsMapData<T> md;
};


} // namespace gismo
