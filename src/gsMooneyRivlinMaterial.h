/** @file gsMooneyRivlinMaterial.h

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

namespace gismo
{

template <class T>
class gsMooneyRivlinMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsMooneyRivlinMaterial( const T c1,
                            const T c2,
                            short_t d)
    :
    gsMooneyRivlinMaterial( gsConstantFunction<T>(c1,d),
                            gsConstantFunction<T>(c2,d))
    {
    }

    gsMooneyRivlinMaterial( const gsFunctionSet<T> & c1,
                            const gsFunctionSet<T> & c2)
    :
    Base()
    {
        this->setParameter(0,c1);
        this->setParameter(1,c2);
    }

    void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const
    {
        GISMO_NO_IMPLEMENTATION;
        // const short_t dim = data.dim;
        // const index_t N = data.size;
        // // Compute the deformation gradient and strain
        // gsMatrix<T> Fres, Eres;
        // Base::eval_deformation_gradient_and_strain_into(patch,u,Fres,Eres);

        // // Resize the result
        // result.resize(dim*dim,u.cols());

        // // Lamé parameters
        // T E, nu;
        // T lambda, mu;
        // gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        // gsMatrix<T> RCG, RCGinv;
        // T J;
        // for (index_t i=0; i!=u.cols(); i++)
        // {
        //     E = m_data.mine().m_parmat(0,i);
        //     nu= m_data.mine().m_parmat(1,i);
        //     lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
        //     mu     = E / ( 2. * ( 1. + nu ) );

        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,dim,dim);
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,dim,dim);
        //     gsAsMatrix<T, Dynamic, Dynamic> S = result.reshapeCol(i,dim,dim);
        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();
        //     S = (lambda*log(J)-mu)*RCGinv + mu*I;
        // }
    }

    void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const
    {
        GISMO_NO_IMPLEMENTATION;
        // const short_t dim = data.dim;
        // const index_t N = data.size;
        // // Compute the deformation gradient and strain
        // gsMatrix<T> Fres, Eres;
        // Base::eval_deformation_gradient_and_strain_into(patch,u,Fres,Eres);

        // // Voigt-size of the tensor
        // const index_t sz = (dim+1)*dim/2;

        // // Resize the result
        // result.resize(sz*sz,u.cols());

        // // Identity tensor
        // gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        // gsMatrix<T> C, Ctemp;
        // gsMatrix<T> RCG, RCGinv;
        // T J;

        // // C tensors
        // gsMatrix<T> Clambda, Cmu;
        // matrixTraceTensor<T>(Clambda,I,I);
        // symmetricIdentityTensor<T>(Cmu,I);

        // // Lamé parameters
        // T E, nu;
        // T lambda, mu;
        // for (index_t i=0; i!=u.cols(); i++)
        // {
        //     E = m_data.mine().m_parmat(0,i);
        //     nu= m_data.mine().m_parmat(1,i);
        //     lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
        //     mu     = E / ( 2. * ( 1. + nu ) );

        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,dim,dim);
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,dim,dim);

        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();

        //     // Compute C
        //     matrixTraceTensor<T>(C,RCGinv,RCGinv);
        //     C *= lambda;
        //     symmetricIdentityTensor<T>(Ctemp,RCGinv);
        //     C += (mu-lambda*log(J))*Ctemp;
        //     result.reshapeCol(i,sz,sz) = C;
        // }
    }
};

}
