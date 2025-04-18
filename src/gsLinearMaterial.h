/** @file gsLinearMaterial.h

    @brief Provides a linear elastic material model
    @todo check equation:
    \f{align*}{
        \Psi(\mathbf{F}) &= \frac{\lambda}{2} \text{tr}^2(\mathbf{F}) + \mu \text{tr}(\mathbf{F}^T \mathbf{F})\\
        \mathbf{S} &= \lambda \text{tr}(\mathbf{E}) \mathbf{I} + 2\mu \mathbf{E}\\
        \mathbf{C} &= \lambda \mathbf{C}^{\text{vol}} + \mu \mathbf{C}^{\text{dev}}
    \f}

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    H.M.Verhelst (2019 - ...., TU Delft)
*/

#pragma once

#include <gsElasticity/src/gsMaterialBase.h>
#include <gsElasticity/src/gsVisitorElUtils.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

/**
 * @brief The gsLinearMaterial class provides a linear elastic material model
 * @ingroup Elasticity
 * @tparam T
 */
template <class T>
class gsLinearMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    /**
     * @brief Constructor with constant parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     * @param dim Dimension of the problem
     */
    gsLinearMaterial(   const T E,
                        const T nu,
                        const size_t dim)
    :
    gsLinearMaterial(gsConstantFunction<T>(E,dim),gsConstantFunction<T>(nu,dim))
    {
    }

    /**
     * @brief Constructor with function parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     */
    gsLinearMaterial(const gsFunctionSet<T> & E,
                     const gsFunctionSet<T> & nu)
    :
    Base()
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    using Base::eval_stress_into;
    /// See \ref gsMaterialBase.h for more details
    void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Resize the result
        Sresult.resize(dim*dim,N);

        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        // Lamé parameters
        gsMatrix<T> lambda(1,N), mu(1,N);
        gsMatrix<T> E     = data.m_pars[0];
        gsMatrix<T> nu    = data.m_pars[1];
        lambda.array() = E.array() * nu.array() / ( ( 1. + nu.array() ) * ( 1. - 2. * nu.array() ) );
        mu.array()     = E.array() / ( 2. * ( 1. + nu.array() ) );
        for (index_t i=0; i!=N; i++)
        {
            gsAsMatrix<T, Dynamic, Dynamic> Emat = data.m_E.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);
            S = lambda(0,i)*Emat.trace()*I + 2*mu(0,i)*Emat;
        }
    }

    using Base::eval_matrix_into;
    /// See \ref gsMaterialBase.h for more details
    void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Voigt-size of the tensor
        const index_t sz = (dim+1)*dim/2;

        // Resize the result
        Cresult.resize(sz*sz,N);

        // Identity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);

        // C tensors
        gsMatrix<T> Clambda, Cmu;
        matrixTraceTensor<T>(Clambda,I,I);
        symmetricIdentityTensor<T>(Cmu,I);

        // Lamé parameters
        gsMatrix<T> lambda(1,N), mu(1,N);
        gsMatrix<T> E     = data.m_pars[0];
        gsMatrix<T> nu    = data.m_pars[1];
        lambda.array() = E.array() * nu.array() / ( ( 1. + nu.array() ) * ( 1. - 2. * nu.array() ) );
        mu.array()     = E.array() / ( 2. * ( 1. + nu.array() ) );
        for (index_t i=0; i!=N; i++)
        {
            // Compute C
            Cresult.reshapeCol(i,sz,sz) = lambda(0,i)*Clambda + mu(0,i)*Cmu;
        }
    }

};

}
