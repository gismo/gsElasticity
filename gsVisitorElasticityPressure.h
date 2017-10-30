/** @file gsVisitorElasticityPressure.h

    @brief Neumann conditions visitor for 2d/3D elasticity.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Shamanskiy
    Inspired by gsVisitorElasticityNeumann by O. Weeger
*/

#pragma once

namespace gismo
{


template <class T>
class gsVisitorElasticityPressure
{
public:

    gsVisitorElasticityPressure(const std::map<unsigned,T> & pressure, boxSide s, gsMatrix<T> & rhsExtra) :
    m_rhsExtra(rhsExtra), pressureData(pressure), side(s)
    {

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
        basis.active_into(quNodes.col(0),localActiveIndices);
        numActiveFunctions = localActiveIndices.rows();
 
        // Evaluate basis functions on element
        basis.eval_into(quNodes, basisValues);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Compute pressure values
        pressureValues.setZero(1,quNodes.cols());
        for (int j = 0; j < numActiveFunctions; ++j)
        {
            if(pressureData.find(localActiveIndices.at(j)) != pressureData.end())
                pressureValues += basisValues.row(j)*pressureData.at(localActiveIndices.at(j));
        }

        // Initialize local matrix/rhs
        localRhs.setZero(m_dim*numActiveFunctions, 1 );
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
            const T weight = quWeights[k] * pressureValues.at(k);

            for (index_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActiveFunctions,numActiveFunctions).noalias() -=
                     weight * unormal(j,0) * basisValues.col(k) ;
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
            gsMatrix<unsigned> ci_actives(localActiveIndices);
            mappers[ci].localToGlobal(localActiveIndices, patchIndex, ci_actives);
			
            for (index_t ai = 0; ai < numActiveFunctions; ++ai)
            {
                const index_t gi = ci * numActiveFunctions +  ai; // row index
                const index_t ii = ci_actives(ai);

                if ( mappers[ci].is_free_index(ii) )
                    m_rhsExtra.row(ii) += localRhs.row(gi);
            }
		}
    }

protected:

	index_t m_dim;
    gsMatrix<T> localRhs; // Local rhs
    gsMatrix <T> & m_rhsExtra; // Global rhs
    
    // Pressure info
    const std::map<unsigned,T> & pressureData; // map of active indices and corresponding values
    boxSide side;
    gsMatrix<T> pressureValues;

    // Basis values
    gsMatrix<T>      basisValues;
    gsMatrix<unsigned> localActiveIndices;
    index_t numActiveFunctions;

    // Normal and Neumann values
	gsVector<T> unormal;
    gsMatrix<T> neuData;   
};


} // namespace gismo
