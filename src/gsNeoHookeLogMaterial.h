/** @file gsNeoHookeLogMaterial.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    O. Weeger    (2012 - 2015, TU Kaiserslautern),
    A.Shamanskiy (2016 - 2020, TU Kaiserslautern),
    H.M.Verhelst (2019 - ...., TU Delft)
*/

#pragma once

#include <gsElasticity/gsMaterialBase.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

template <class T>
class gsNeoHookeLogMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsNeoHookeLogMaterial(  const T E,
                            const T nu,
                            short_t d)
    :
    gsNeoHookeLogMaterial(gsConstantFunction<T>(E,d),
                     gsConstantFunction<T>(nu,d))
    {
    }

    gsNeoHookeLogMaterial(  const gsFunctionSet<T> & E,
                            const gsFunctionSet<T> & nu)
    :
    Base()
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    gsNeoHookeLogMaterial(  const T E,
                            const T nu,
                            const gsFunctionSet<T> & patches)
    :
    gsNeoHookeLogMaterial(gsConstantFunction<T>(E,patches.domainDim()),
                          gsConstantFunction<T>(nu,patches.domainDim()),patches)
    {
    }

    gsNeoHookeLogMaterial(  const gsFunctionSet<T> & E,
                            const gsFunctionSet<T> & nu,
                            const gsFunctionSet<T> & patches)
    :
    Base(&patches)
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const override
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
            E = data.m_parmat(0,i);
            nu= data.m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> F = data.m_F.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);
            J = F.determinant();
            RCG = F.transpose() * F;
            RCGinv = RCG.cramerInverse();
            S = (lambda*log(J)-mu)*RCGinv + mu*I;
        }
    }

    void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const override
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
            E = data.m_parmat(0,i);
            nu= data.m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> F = data.m_F.reshapeCol(i,dim,dim);

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
