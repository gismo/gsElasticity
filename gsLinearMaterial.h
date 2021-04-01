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

namespace gismo
{

template <class T>
class gsLinearMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsLinearMaterial(   const T E,
                        const T nu)
    {
        m_Emodulus = E;
        m_PoissonRatio = nu;
    }

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        const index_t sz = (u.rows()+1)*u.rows()/2;
        result.resize(sz*sz,u.cols());

        gsMatrix<T> C, Ctemp;
        T lambda = m_Emodulus * m_PoissonRatio / ( ( 1. + m_PoissonRatio ) * ( 1. - 2. * m_PoissonRatio ) );
        T mu     = m_Emodulus / ( 2. * ( 1. + m_PoissonRatio ) );
        // linear elasticity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
        matrixTraceTensor<T>(C,I,I);
        C *= lambda;
        symmetricIdentityTensor<T>(Ctemp,I);
        C += mu*Ctemp;

        for (index_t k=0; k!=u.cols(); k++)
            result.reshapeCol(k,sz,sz) = C;

    }

    /// \a u points, \a k patch index
    void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    T m_Emodulus, m_PoissonRatio;

};

}
