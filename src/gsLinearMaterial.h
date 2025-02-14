/** @file gsLinearMaterial.h

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
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsUtils/gsThreaded.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

template <class T>
class gsLinearMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsLinearMaterial(   const T E,
                        const T nu,
                        const gsFunctionSet<T> & patches,
                        const gsFunctionSet<T> & deformed)
    :
    gsLinearMaterial(gsConstantFunction<T>(E,patches.domainDim()),
                     gsConstantFunction<T>(nu,patches.domainDim()),patches,deformed)
    {
    }

    gsLinearMaterial(const gsFunctionSet<T> & E,
                     const gsFunctionSet<T> & nu,
                     const gsFunctionSet<T> & patches,
                     const gsFunctionSet<T> & deformed)
    :
    Base(&patches,&deformed)
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    gsLinearMaterial(   const T E,
                        const T nu,
                        const gsFunctionSet<T> & patches)
    :
    gsLinearMaterial(gsConstantFunction<T>(E,patches.domainDim()),
                     gsConstantFunction<T>(nu,patches.domainDim()),patches)
    {
    }

    gsLinearMaterial(const gsFunctionSet<T> & E,
                     const gsFunctionSet<T> & nu,
                     const gsFunctionSet<T> & patches)
    :
    Base(&patches)
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    gsLinearMaterial(   const T E,
                        const T nu,
                        const size_t dim)
    :
    gsLinearMaterial(gsConstantFunction<T>(E,dim),gsConstantFunction<T>(nu,dim))
    {
    }

    gsLinearMaterial(const gsFunctionSet<T> & E,
                     const gsFunctionSet<T> & nu)
    :
    Base()
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Resize the result
        Sresult.resize(dim*dim,N);

        // Lamé parameters
        T E, nu;
        T lambda, mu;
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);

        for (index_t i=0; i!=N; i++)
        {
            E = m_data.mine().m_parmat(0,i);
            nu= m_data.mine().m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> E = data.m_E.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);
            S = lambda*E.trace()*I + 2*mu*E;
        }
    }

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
        T E, nu;
        T lambda, mu;
        for (index_t i=0; i!=N; i++)
        {
            E = m_data.mine().m_parmat(0,i);
            nu= m_data.mine().m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );
            // Compute C
            Cresult.reshapeCol(i,sz,sz) = lambda*Clambda + mu*Cmu;
        }
    }

protected:
    using Base::m_data;

};

}
