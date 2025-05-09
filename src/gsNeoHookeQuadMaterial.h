/** @file gsNeoHookeQuadMaterial.h

    @brief Provides a neo-Hookean material model
    @todo equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    H.M.Verhelst (2019 - ...., TU Delft)
*/

#pragma once

#include <gsElasticity/src/gsMaterialBase.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

/**
 * @brief The gsNeoHookeQuadMaterial class provides a neo-Hookean material model
 * @ingroup Elasticity
 * @tparam T
 */
template <class T>
class gsNeoHookeQuadMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    /**
     * @brief Constructor with constant parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     * @param dim Dimension of the problem
     */
    gsNeoHookeQuadMaterial( const T E,
                            const T nu,
                            short_t d)
    :
    gsNeoHookeQuadMaterial(gsConstantFunction<T>(E,d),
                     gsConstantFunction<T>(nu,d))
    {
    }

    /**
     * @brief Constructor with function parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     */
    gsNeoHookeQuadMaterial( const gsFunctionSet<T> & E,
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

        // Lamé parameters
        gsMatrix<T> lambda(1,N), mu(1,N);
        gsMatrix<T> E     = data.m_pars[0];
        gsMatrix<T> nu    = data.m_pars[1];
        lambda.array() = E.array() * nu.array() / ( ( 1. + nu.array() ) * ( 1. - 2. * nu.array() ) );
        mu.array()     = E.array() / ( 2. * ( 1. + nu.array() ) );

        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        gsMatrix<T> RCG, RCGinv;
        T J;
        for (index_t i=0; i!=N; i++)
        {
            gsAsMatrix<T, Dynamic, Dynamic> F = data.m_F.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);
            J = F.determinant();
            RCG = F.transpose() * F;
            RCGinv = RCG.cramerInverse();
            S = ((lambda(0,i))*(J*J-1)/2-mu(0,i))*RCGinv + mu(0,i)*I;
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
        gsMatrix<T> C, Ctemp;
        gsMatrix<T> RCG, RCGinv;
        T J;

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
            gsAsMatrix<T, Dynamic, Dynamic> F = data.m_F.reshapeCol(i,dim,dim);

            J = F.determinant();
            RCG = F.transpose() * F;
            RCGinv = RCG.cramerInverse();

            // Compute C
            matrixTraceTensor<T>(C,RCGinv,RCGinv);
            C *= lambda(0,i)*J*J;
            symmetricIdentityTensor<T>(Ctemp,RCGinv);
            C += (mu(0,i)-lambda(0,i)*(J*J-1)/2)*Ctemp;
            Cresult.reshapeCol(i,sz,sz) = C;
        }
    }
};

}
