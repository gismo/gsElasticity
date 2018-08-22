/** @file gsVisitorElThermoBoundary.h

    @brief Computes boundary thermal stress in the normal direction.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/

#pragma once

namespace gismo
{

template <class T>
class gsVisitorElThermoBoundary
{
public:

    gsVisitorElThermoBoundary(const gsFunction<T> & heatSol, boxSide s, gsMatrix<T> & rhsExtra, T initTemp,
                              T lambda, T mu, T thExpCoef) :
        thermoSol(heatSol),side(s),m_rhsExtra(rhsExtra),startTemp(initTemp),
        m_lambda(lambda),m_mu(mu),m_thExpCoef(thExpCoef)
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
        evFlags = NEED_VALUE | NEED_JACOBIAN; // Jacobian?
    }

    inline void evaluate(gsBasis<T> const       & basis,
                         gsGeometryEvaluator<T> & geoEval,
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

        // Compute heat(thermo) values
        thermoSol.eval_into(quNodes,thermoValues);

        // Initialize local rhs
        localRhs.setZero(m_dim*numActiveFunctions, 1 );
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);

            // Multiply quadrature weight by the measure of normal and time-dependent factor
            const T weight = m_thExpCoef*(2*m_mu+m_dim*m_lambda)*quWeights[k] * (thermoValues.at(k)-startTemp);

            for (index_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActiveFunctions,numActiveFunctions).noalias() +=
                     weight * unormal(j,0) * basisValues.col(k) ;
        }
    }

    void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int             patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        for (index_t d = 0; d!= m_dim; ++d)
        {
            gsMatrix<unsigned> globalActiveIndices(localActiveIndices);
            mappers[d].localToGlobal(localActiveIndices, patchIndex, globalActiveIndices);

            for (index_t i = 0; i < numActiveFunctions; ++i)
            {
                const index_t gi = d * numActiveFunctions +  i; // row index
                const index_t ii = globalActiveIndices(i);

                if ( mappers[d].is_free_index(ii) )
                    m_rhsExtra.row(ii) += localRhs.row(gi);
            }
        }
    }

protected:

    // Thermo info
    const gsFunction<T> & thermoSol;
    boxSide side;
    gsMatrix<T> thermoValues;

    index_t m_dim;
    gsMatrix<T> localRhs; // Local rhs
    gsMatrix <T> & m_rhsExtra; // Global rhs

    T startTemp;
    T m_lambda;
    T m_mu;
    T m_thExpCoef;


    // Basis values
    gsMatrix<T>      basisValues;
    gsMatrix<unsigned> localActiveIndices;
    index_t numActiveFunctions;

    // Normal values
    gsVector<T> unormal;



}; //class definition ends

} // namespace ends
