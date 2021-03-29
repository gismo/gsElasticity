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
                        const T nu,
                        const index_t dim
                                        )
    {
        m_Emodulus = E;
        m_PoissonRatio = nu;
        Base::m_dim = dim;
        if (dim==2)
            Base::m_size = 3;
        else if (dim==3)
            Base::m_size = 6;
        else
            GISMO_ERROR("Material matrix for dimension "<<dim<<" unknown");
    }

    /// Implementation of eval_into, see \ref gsFunction
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(Base::m_size*Base::m_size,u.cols());

        gsMatrix<T> C, Ctemp;
        T lambda = m_Emodulus * m_PoissonRatio / ( ( 1. + m_PoissonRatio ) * ( 1. - 2. * m_PoissonRatio ) );
        T mu     = m_Emodulus / ( 2. * ( 1. + m_PoissonRatio ) );
        // linear elasticity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(Base::m_dim,Base::m_dim);
        matrixTraceTensor<T>(C,I,I);
        C *= lambda;
        symmetricIdentityTensor<T>(Ctemp,I);
        C += mu*Ctemp;

        for (index_t k=0; k!=u.cols(); k++)
            result.reshapeCol(k,Base::m_size,Base::m_size) = C;

    }

    /// Implementation of eval_into, see \ref gsFunction
    void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    T m_Emodulus, m_PoissonRatio;

};

}