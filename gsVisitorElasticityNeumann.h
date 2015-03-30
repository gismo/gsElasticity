/** @file gsVisitorElasticityNeumann.h

    @brief Neumann conditions visitor for 2d/3D elasticity.

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
class gsVisitorElasticityNeumann
{
public:

    gsVisitorElasticityNeumann(const gsFunction<T> & neudata, boxSide s, T tfac = 1.0) : 
    neudata_ptr(&neudata), side(s)
    { 
        m_tfac = tfac;
	}

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        m_dim = basis.dim();
		
		const int dir = side.direction();
        gsVector<int> numQuadNodes ( m_dim );
        for (int i = 0; i < m_dim; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         // todo: add element here for efficiency
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(quNodes.col(0) , actives);
        numActive = actives.rows();
 
        // Evaluate basis functions on element
        basis.eval_into(quNodes, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Neumann data
        neudata_ptr->eval_into(geoEval.values(), neuData);

        // Initialize local matrix/rhs
        localRhs.setZero(m_dim*numActive, 1 );
    }

    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            /*
			// Compute the outer normal vector on the side
            //geoEval.outerNormal(k, side, unormal); // (!)
            //geoEval.normal(k,unormal);  
            
            // Multiply quadrature weight by the measure of normal on the boundary
            const T weight = quWeights[k] * geoEval.jacobian(k).col( ! side.direction() ).norm();
            */

			// Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);
            
            // Multiply quadrature weight by the measure of normal and time-dependent factor
            const T weight = quWeights[k] * unormal.norm() * m_tfac;

            for (index_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                    weight * neuData(j,k) * basisData.col(k) ;
        }
    }
    
    void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int             patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
		for (index_t ci = 0; ci!= m_dim; ++ci)
		{          
			gsMatrix<unsigned> ci_actives(actives);
			mappers[ci].localToGlobal(actives, patchIndex, ci_actives);
			
			for (index_t ai = 0; ai < numActive; ++ai)
            {
                const index_t gi = ci * numActive +  ai; // row index
                const index_t ii = ci_actives(ai);

                if ( mappers[ci].is_free_index(ii) )
                    rhsMatrix.row(ii) += localRhs.row(gi);
            }
		}
    }

protected:

    index_t numActive;
	index_t m_dim;
    
    // Neumann function
    const gsFunction<T> * neudata_ptr;
    boxSide side;
	T m_tfac;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<unsigned> actives;

    // Normal and Neumann values
	gsVector<T> unormal;
    gsMatrix<T> neuData;

    // Local matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
