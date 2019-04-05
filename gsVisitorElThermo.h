/** @file gsVisitorElThermo.h

    @brief Computes thermal stress in the interia.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/
/*
#pragma once

namespace gismo
{

template <class T>
class gsVisitorElThermo
{
public:

    gsVisitorElThermo(const gsFunction<T> & heatSol, gsMatrix<T> & rhsExtra,
                      T lambda, T mu, T thExpCoef) :
        thermoSol(heatSol),m_rhsExtra(rhsExtra),
        m_lambda(lambda),m_mu(mu),m_thExpCoef(thExpCoef)
    {

    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        m_dim = basis.dim();
        gsVector<index_t> numQuadNodes(m_dim);

        for (index_t i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

    }

    inline void evaluate(gsBasis<T> const       & basis,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        basis.active_into(quNodes.col(0), localActiveIndices);
        numActiveFunctions = localActiveIndices.rows();

        // Evaluate basis functions on element
        basis.evalAllDers_into(quNodes, 1, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);

        heatGrad.setZero(2,quNodes.cols());
        localRhs.setZero(m_dim*numActiveFunctions,1);
        thermoSol.deriv_into(quNodes,heatGrad);
    }

    inline void assemble(gsDomainIterator<T>    & , // element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k,heatGrad,physGrad);

            const T weight = m_thExpCoef*(2*m_mu+m_dim*m_lambda)*quWeights[k]*geoEval.measure(k);

            for (index_t d = 0; d < m_dim; ++d)
            {
                localRhs.middleRows(d*numActiveFunctions,numActiveFunctions).noalias() -=
                        weight * physGrad(d,0) * basisData[0].col(k);
            }
        }
    }

    void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                       const gsMatrix<T>     & , // eliminatedDofs,
                       const int               patchIndex,
                       gsSparseMatrix<T>     & , // sysMatrix,
                       gsMatrix<T>           & ) // rhsMatrix )
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

    const gsFunction<T> & thermoSol;
    gsMatrix<T> heatGrad, physGrad;

    index_t m_dim;
    gsMatrix<T> localRhs; // Local rhs
    gsMatrix <T> & m_rhsExtra; // Global rhs

    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<unsigned> localActiveIndices;
    index_t numActiveFunctions;

    T m_lambda;
    T m_mu;
    T m_thExpCoef;



}; //class definition ends

} // namespace ends
*/
