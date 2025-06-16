/** @file gsNeoHookeLogMaterial.h

    @brief Provides a neo-Hookean material model with logarithmic strain
    @todo check equation:
    \f{align*}{
        \Psi(\mathbf{F}) &= \frac{\lambda}{2} \log^2(J) - \mu \log(J) + \frac{\mu}{2} \text{tr}(\mathbf{C}^{\text{vol}})\\
        \mathbf{S} &= \lambda \log(J) \mathbf{C}^{-1} + \mu \mathbf{I}\\
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

#include <gsElasticity/gsMaterialBase.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

/**
 * @brief The gsNeoHookeLogMaterial class provides a neo-Hookean material model with logarithmic strain
 * @ingroup Elasticity
 * @tparam T
 */
template <class T>
class gsNeoHookeLogMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    /**
     * @brief Constructor with constant parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     * @param dim Dimension of the problem
     */
    gsNeoHookeLogMaterial(  const T E,
                            const T nu,
                            short_t d)
    :
    gsNeoHookeLogMaterial(gsConstantFunction<T>(E,d),
                     gsConstantFunction<T>(nu,d))
    {
    }

    /**
     * @brief Constructor with function parameters
     * @param E Young's modulus
     * @param nu Poisson's ratio
     */
    gsNeoHookeLogMaterial(  const gsFunctionSet<T> & E,
                            const gsFunctionSet<T> & nu)
    :
    Base()
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    /// See \ref gsMaterialBase.h for more details
    void compute_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const override
    {
        this->eval_stress_into(data, Sresult);
    }

    /// See \ref gsMaterialBase.h for more details
    static void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult)
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Resize the result
        Sresult.resize(dim*dim,N);

        // Lamé parameters
        T E, nu;
        T lambda, mu;
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        gsMatrix<T> RCG, RCGinv;
        T J;
        for (index_t i=0; i!=N; i++)
        {
            E = data.parameters[0](0,i);
            nu= data.parameters[1](0,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> F = data.deformationGradient.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);
            J = F.determinant();
            RCG = F.transpose() * F;
            RCGinv = RCG.cramerInverse();
            S = (lambda*log(J)-mu)*RCGinv + mu*I;
        }
    }


    /// See \ref gsMaterialBase.h for more details
    void compute_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const override
    {
        this->eval_matrix_into(data, Sresult);
    }

    /// See \ref gsMaterialBase.h for more details
    static void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult)
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
        T E, nu;
        T lambda, mu;
        for (index_t i=0; i!=N; i++)
        {
            E = data.parameters[0](0,i);
            nu= data.parameters[1](0,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> F = data.deformationGradient.reshapeCol(i,dim,dim);

            J = F.determinant();
            RCG = F.transpose() * F;
            RCGinv = RCG.cramerInverse();

            // Compute C
            matrixTraceTensor<T>(C,RCGinv,RCGinv);
            C *= lambda;
            symmetricIdentityTensor<T>(Ctemp,RCGinv);
            C += (mu-lambda*log(J))*Ctemp;
            Cresult.reshapeCol(i,sz,sz) = C;
        }
    }
};

}
