/** @file gsResidualLinElast.h

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

/** \brief Computes residual for the linear elasticity problem.
 *
 * The problem is given by
 * \f[ - \mathrm{ div }\ \sigma(u) = f \mathrm{~on~}\Omega,\quad u = u_D \mathrm{~on~} \partial \Omega, \f]
 * where \n
 * \f$ \sigma \f$ is the Cauchy-stress-tensor\n
 * \f$ f \f$ are given body forces\n
 * \f$ u_D \f$ are given displacement boundary conditions\n
 * \f$ \Omega \f$ is the computational domain
 *
 * The residual is comuted as
 * \f[ \| \mathrm{ div }\ \sigma(u) + f \|_{L_2} \f]
 *
 * \remarks Assumes that the dimensions of the parameter space
 * and the target space are the same
 * (i.e., \f$ u: \mathbb R^d \rightarrow \mathbb R^d\f$).
 *
 * \warning Note that the terms regarding Neumann (traction) boundary
 * conditions
 * and jumps of the normal derivative across element interfaces are
 * NEGLECTED in the current version (13.Jul.2015).
 *
 * \warning Work in progress!
 *
 * \ingroup Assembler
 */
template <class T>
class gsResidualLinElast : public gsNorm<T>
{
    friend class gsNorm<T>;

public:

    /** \brief Constructor
     * \param _discSolution Discrete solution (displacement field)
     * \param E_modulus Young's modulus.
     * \param poissons_ratio Poisson's ratio.
     * \param _bodyForces Right-hand-side-/body-forces-function \f$ f\f$.
     * \param _rhsFunctionParam Flag indicating whether the \em _rhsFunction
     * is parameterized (in this case, the evaluation points must be given
     * on the parameter domain
     */
    gsResidualLinElast(const gsField<T> & _discSolution,
                       const T E_modulus,
                       const T poissons_ratio,
                       const gsFunction<T> & _bodyForces,
                       bool _exactSolutionParam = false)
    : gsNorm<T>(_discSolution,_bodyForces), m_exSolParam(_exactSolutionParam)
    {
        m_lambda = E_modulus * poissons_ratio / ( (1.+poissons_ratio)*(1.-2.*poissons_ratio)) ;
        m_mu     = E_modulus / (2.*(1.+poissons_ratio)) ;
    }

public:

    /** \brief Computes the error estimate.
     *
     * Computes the residual-based error estimate \f$\eta\f$
     * (see class-documentation at the top).
     *
     * \param storeEltWise Bool indicating whether the element-wise
     * errors should be stored also. If <em>storeEletWise = true</em>,
     * the gsVector of element-wise estimates \f$\eta_K\f$
     * can be obtained by
     * calling elementNorms().
     *
     * \returns The total estimated error \f$ \eta \f$.
     *
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
        m_dim2nds = (m_parDim*(m_parDim+1))/2;

        GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Dimension must be 2 or 3.");

        // Setup Quadrature
        gsVector<index_t> numQuadNodes( m_parDim );
        for (unsigned i = 0; i < m_parDim; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN| NEED_2ND_DER| NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsGeometry<T>    & discSolution,
                         const gsFunction<T>    & exactGradients,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate discrete solution
        discSolution.deriv_into(quNodes, m_discSolGrads);
        discSolution.deriv2_into(quNodes, m_discSol2nd);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand-side function (defined of physical domain)
        exactGradients.eval_into(geoEval.values(), m_exactBodyForces);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);

        for (index_t qk = 0; qk < quWeights.size(); ++qk) // loop over quadrature nodes
        {
            const T weight = quWeights[qk] * geoEval.measure(qk);

            //const typename gsMatrix<T>::constColumns J = geoEval.jacobian(qk);
            //gsMatrix<T> Jinv = J.inverse();
            gsMatrix<T> disc2nd;

            // TODO: CHECK IF THIS REALLY IS THE CORRECT FUNCTION TO CALL.
            geoEval.transformDeriv2Hgrad( qk, m_discSolGrads, m_discSol2nd, disc2nd);

            std::vector< gsMatrix<T> > du;
            gsMatrix<T> tdu( m_parDim, m_parDim );
            for( unsigned i = 0; i < m_parDim; i++)
            {
                tdu.setZero();
                if( m_parDim == 2 )
                {
                    tdu(0,0) = disc2nd( i, 0 );
                    tdu(1,1) = disc2nd( i, 1 );
                    tdu(0,1) = disc2nd( i, 2 );
                    tdu(1,0) = disc2nd( i, 2 );
                }
                else if( m_parDim == 3 )
                {
                    tdu(0,0) = disc2nd( i, 0 );
                    tdu(1,1) = disc2nd( i, 1 );
                    tdu(2,2) = disc2nd( i, 2 );
                    tdu(0,1) = disc2nd( i, 3 );
                    tdu(1,0) = disc2nd( i, 3 );
                    tdu(0,2) = disc2nd( i, 4 );
                    tdu(2,0) = disc2nd( i, 4 );
                    tdu(1,2) = disc2nd( i, 5 );
                    tdu(2,1) = disc2nd( i, 5 );
                }
                du.push_back( tdu );
            }

            T diff(0.0);
            for( unsigned i=0; i < m_parDim; i++)
            {
                T Li = T(0.0);
                for( unsigned j=0; j < m_parDim; j++)
                    Li += m_lambda * du[j](j,i) + m_mu * ( du[j](i,j) + du[i](j,j) );
                diff += ( Li - m_exactBodyForces(i,qk) )*( Li - m_exactBodyForces(i,qk) );
            }

            sum += weight * diff;
        } // qk

        return sum;
    }

private:

    gsMatrix<T> m_discSolGrads, m_discSol2nd, m_exactBodyForces;
    gsMatrix<T> m_C;
    T m_mu;
    T m_lambda;

    size_t m_dim2nds;

    unsigned m_parDim;

    bool m_exSolParam;
}; // class


} // namespace gismo

