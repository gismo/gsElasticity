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

#include <gsElasticity/gsMaterialBaseDim.h>
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

template <short_t d, class T>
class gsLinearMaterial : public gsMaterialBaseDim<d,T>
{

public:
    using Base = gsMaterialBaseDim<d,T>;

    gsLinearMaterial(   const T E,
                        const T nu,
                        const gsFunctionSet<T> & patches,
                        const gsFunctionSet<T> & deformed)
    :
    gsLinearMaterial(gsConstantFunction<T>(E,d),gsConstantFunction<T>(nu,d),patches,deformed)
    {
    }

    gsLinearMaterial(const gsFunctionSet<T> & E,
                     const gsFunctionSet<T> & nu,
                     const gsFunctionSet<T> & patches,
                     const gsFunctionSet<T> & deformed)
    :
    Base(&patches,&deformed,nullptr)
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
    }

    void eval_stress_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u);
        // Compute the strain
        gsMatrix<T> Eres;
        Base::eval_strain_into(patch,u,Eres);

        // Lamé parameters
        T E, nu;
        T lambda, mu;
        gsMatrix<T,d,d> I = gsMatrix<T>::Identity(d,d);
        for (index_t i=0; i!=u.cols(); i++)
        {
            E = m_data.mine().m_parmat(0,i);
            nu= m_data.mine().m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );

            gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,d,d);
            result.reshapeCol(i,d,d) = lambda*E.trace()*I + 2*mu*E;
        }

    }

    void eval_matrix_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u);
        // Compute the strain
        gsMatrix<T> Eres;
        Base::eval_strain_into(patch,u,Eres);

        // Voight-size of the tensor
        const index_t sz = (d+1)*d/2;

        // Resize the result
        result.resize(sz*sz,u.cols());

        // Identity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(d,d);

        // C tensors
        gsMatrix<T> Clambda, Cmu;
        matrixTraceTensor<T>(Clambda,I,I);
        symmetricIdentityTensor<T>(Cmu,I);

        // Lamé parameters
        T E, nu;
        T lambda, mu;
        for (index_t i=0; i!=u.cols(); i++)
        {
            E = m_data.mine().m_parmat(0,i);
            nu= m_data.mine().m_parmat(1,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );
            // Compute C
            result.reshapeCol(i,sz,sz) = lambda*Clambda + mu*Cmu;
        }
    }

protected:
    using Base::m_data;

};

}
