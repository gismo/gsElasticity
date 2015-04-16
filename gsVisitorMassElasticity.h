/** @file gsVisitorMassElasticity.h

    @brief Element visitor for elasticity mass matrix for 2D plain strain and 3D continua.

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
class gsVisitorMassElasticity
{
public:

    /// Constructor
    gsVisitorMassElasticity(T rho) : 
	m_rho(rho)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        m_dim = basis.dim();

		gsVector<index_t> numQuadNodes(m_dim);
        for (size_t i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE;
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
        basis.eval_into(quNodes, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {           
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);

			for (index_t i = 0; i < numActive; i++)
			{
				// Exploit symmetry and local diagonal form of M
				for (index_t j = i; j < numActive; j++)
				{
					const T int_val = weight * m_rho * basisData(i,k) * basisData(j,k);
						
					localMat(i, j) += int_val;
				}
			}

        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }
    
    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const int patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
						      gsMatrix<T>           & rhsMatrix
							  )
    {
       	
		// Local DoFs to global DoFs
		std::vector< gsMatrix<unsigned> > ci_actives(m_dim,actives);

        for (size_t ci = 0; ci != m_dim; ++ci)
			mappers[ci].localToGlobal(actives, patchIndex, ci_actives[ci]);

        for (index_t ai=0; ai < numActive; ++ai)    
		{
			const index_t gi = ai; // row index

			for (size_t ci = 0; ci!= m_dim; ++ci)   
            {               
                const index_t ii = ci_actives[ci](ai);

                if ( mappers[ci].is_free_index(ii) )
                {                   
					// Exploit symmetry of M 
					for (index_t aj=ai; aj < numActive; ++aj)
					{
                        const index_t gj = aj; // column index

						// Exploit local diagonal form of M 
						size_t cj = ci;
                       
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
    
    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    // Basis values
    gsMatrix<T>        basisData;
    gsMatrix<unsigned> actives;
    index_t            numActive;

protected:

	// Dimension
	size_t m_dim;

    // Rho
    T m_rho;

    // Local matrices
    gsMatrix<T> localMat;
};


} // namespace gismo

