/** @file gsEnergyNormLinElast.h

    @brief Energy-norm for the linear elasticity problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#pragma once

#include <gsAssembler/gsNorm.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** \brief Computes the error \f$ u_h - v \f$ in the energy norm
 *  for linear elasticity problems
 *
 * ...which is given by
 * \n
 * \f[ \| u_h - v \|_E := \sqrt{ a(u_h-v, u_h-v )}, \f]
 * where\n
 * \f$ a(\cdot,\cdot) \f$ is the bilinear form induced by
 * the linear elasticity problem\n
 * \f$ u_h \f$ is a computed discrete solution\n
 * \f$ v \f$ is a given function (e.g., a known exact solution)\n
 *
 * By setting \f$ v = 0\f$, one can compute the energy norm
 * of \f$ u_h \f$.
 *
 * \remarks Assumes that the dimensions of the parameter space
 * and the target space are the same
 * (i.e., \f$ u_h: \mathbb R^d \rightarrow \mathbb R^d\f$).
 *
 * \warning Work in Progress!
 *
 * \ingroup Assembler
 */
template <class T>
class gsEnergyNormLinElast : public gsNorm<T>
{
    friend class gsNorm<T>;

public:

    /** \brief Constructor
     *
     * \param _discSolution Discrete solution (displacement field)
     * \param E_modulus Young's modulus.
     * \param poissons_ratio Poisson's ratio.
     * \param _exactGradients A gsFunction providing the exact \b gradients
     * of the function v.
     *
     * \param _exactGradsParam Flag indicating whether the \em _exactGradients
     * are parameterized (in this case, the evaluation points must be given
     * on the parameter domain). Defaults to \em false, which means that
     * \em _exactGradients can be evaluated on the physical domain.
     */
    gsEnergyNormLinElast(const gsField<T> & _discSolution,
                         const T E_modulus,
                         const T poissons_ratio,
                         const gsFunction<T> & _exactGradients,
                         bool _exactGradsParam = false)
    : gsNorm<T>(_discSolution,_exactGradients), m_exSolParam(_exactGradsParam)
    {
        m_lambda = E_modulus * poissons_ratio / ( (1.+poissons_ratio)*(1.-2.*poissons_ratio)) ;
        m_mu     = E_modulus / (2.*(1.+poissons_ratio)) ;
    }

public:

    /** \brief Computes the error.
     *
     * Computes the error in the energy norm
     * (see class-documentation at the top).
     *
     *
     * \param storeEltWise Bool indicating whether the element-wise
     * errors should be stored also. If <em>storeEletWise = true</em>,
     * the gsVector of element-wise estimates \f$\eta_K\f$
     * can be obtained by
     * calling elementNorms().
     *
     * \returns The total error \f$ \eta \f$.
     */
    T compute(bool storeElWise = true)
    {
        this->apply(*this,storeElWise);
        return this->m_value;
    }


protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        m_parDim = basis.dim();
        m_dimStr = (m_parDim*(m_parDim+1))/2;

        GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Dimension must be 2 or 3.");

        // Setup Quadrature
        gsVector<index_t> numQuadNodes( m_parDim );
        for (unsigned i = 0; i < m_parDim; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsGeometry<T>    & discSolution,
                         const gsFunction<T>    & exactGradients,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate discrete solution
        discSolution.deriv_into(quNodes, m_discSolGrads);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand-side function (defined of physical domain)
        exactGradients.eval_into(geoEval.values(), m_exactSolGrads);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);

        for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);

            const typename gsMatrix<T>::constColumns J = geoEval.jacobian(k);
            gsMatrix<T> Jinv = J.inverse();
            gsMatrix<T> discSolG;

            discSolG = m_discSolGrads.col(k);
            discSolG.resize( m_parDim, m_parDim );
            discSolG = discSolG * Jinv.transpose();

            gsMatrix<T> exSolG = m_exactSolGrads.col(k);
            exSolG.resize(m_parDim, m_parDim );

            gsMatrix<T> useGr = discSolG - exSolG;

            gsVector<T> sT( m_dimStr );

            // the following could probably done with more memory-efficiency
            // but for the debugging-process, readability had higher priority.
            T tmpNorm(0.0);
            if( m_parDim == 2 )
            {
                sT[0] = useGr(0,0);
                sT[1] = useGr(1,1);
                sT[2] = ( useGr(0,1) + useGr(1,0) );

                T tmpTr = m_lambda*(sT[0]+sT[1]);
                tmpNorm = ( tmpTr + 2*m_mu*sT[0] ) * sT[0]
                        + ( tmpTr + 2*m_mu*sT[1] ) * sT[1]
                        + m_mu * sT[2] * sT[2];

            }
            else
            {
                sT[0] = useGr(0,0);
                sT[1] = useGr(1,1);
                sT[2] = useGr(2,2);
                sT[3] = ( useGr(0,1) + useGr(1,0) );
                sT[4] = ( useGr(0,2) + useGr(2,0) );
                sT[5] = ( useGr(1,2) + useGr(2,1) );

                T tmpTr = m_lambda*(sT[0]+sT[1]+sT[2]);

                tmpNorm = ( tmpTr + 2*m_mu*sT[0] ) * sT[0]
                        + ( tmpTr + 2*m_mu*sT[1] ) * sT[1]
                        + ( tmpTr + 2*m_mu*sT[2] ) * sT[2]
                        + m_mu * sT[3] * sT[3]
                        + m_mu * sT[4] * sT[4]
                        + m_mu * sT[5] * sT[5];
            }
            sum += weight * tmpNorm;
        }
        return sum;
    }

private:

    gsMatrix<T> m_discSolGrads, m_discSolGradsPhys, m_exactSolGrads;
    gsMatrix<T> m_C;
    T m_mu;
    T m_lambda;

    size_t m_dimStr;
    //size_t m_dim;

    unsigned m_parDim;

    bool m_exSolParam;
}; // class


} // namespace gismo

