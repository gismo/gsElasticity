/** @file gsLinearDegradedMaterial.h

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
class gsLinearDegradedMaterial : public gsMaterialBaseDim<d,T>
{

public:
    using Base = gsMaterialBaseDim<d,T>;

    gsLinearDegradedMaterial(   const T E,
                                const T nu,
                                const gsFunctionSet<T> & patches,
                                const gsFunctionSet<T> & deformed,
                                const gsFunctionSet<T> & damage)
    :
    gsLinearDegradedMaterial(gsConstantFunction<T>(E,d),gsConstantFunction<T>(nu,d),patches,deformed,damage) //?? anadi damage?
    {
    }

    gsLinearDegradedMaterial(   const gsFunctionSet<T> & E,
                                const gsFunctionSet<T> & nu,
                                const gsFunctionSet<T> & patches,
                                const gsFunctionSet<T> & deformed,
                                const gsFunctionSet<T> & damage)
    :
    Base(&patches,&deformed,nullptr)    
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
        this->setParameter(2,damage);
    }

    void eval_stress_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u); // compute functions (E,u,damage)
        // Compute the strain
        gsMatrix<T> Eres;
        Base::eval_strain_into(patch,u,Eres);

        // Lamé parameters
        T E, nu;
        T lambda, mu, kappa;
        // Damage variable ("parameter")
        T damage;

        gsMatrix<T,d,d> I = gsMatrix<T>::Identity(d,d);
        result.resize(d*d,u.cols());

        for (index_t i=0; i!=u.cols(); i++)
        {
            E      = m_data.mine().m_parmat(0,i);
            nu     = m_data.mine().m_parmat(1,i);
            damage = m_data.mine().m_parmat(2,i);

            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );
            kappa  = E / ( 3. * ( 1. - 2.*nu ) );

            // Compute strain Heaviside function
            gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,d,d);

            auto E_vol = E.trace(); //volumetric strain (scalar) 
            auto E_vol_pos = (E_vol > 0) ? E_vol : 0; // positive Heaviside function
            auto E_vol_neg = (E_vol < 0) ? E_vol : 0; // negative Heaviside function
            auto E_dev = E - 1./d * E_vol * I; // deviatoric strain
            
            // Compute stress (Tensor dxd)
            result = math::pow((1. - damage),2) * (kappa*E_vol_pos*I + 2.*mu*E_dev) + kappa*E_vol_neg*I;
        }

    }

    void eval_matrix_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u);
        
        gsMatrix<T> Eres;
        Base::eval_strain_into(patch,u,Eres);

        // Voigt-size of the tensor
        const index_t sz = (d+1)*d/2;

        // Resize the result
        result.resize(sz*sz,u.cols());

        // Identity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(d,d);

        // C tensors
        gsMatrix<T> Cvol, Cdev;
        matrixTraceTensor<T>(Cvol,I,I);
        deviatoricTensor<T>(Cdev,I);

        // Lamé parameters
        T E, nu;
        T lambda, mu, kappa;
        // Damage parameter
        T damage;

        bool H_pos, H_neg;
        T E_vol;

        for (index_t i=0; i!=u.cols(); i++) 
        {
            E      = m_data.mine().m_parmat(0,i);
            nu     = m_data.mine().m_parmat(1,i);
            damage = m_data.mine().m_parmat(2,i);
            
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            mu     = E / ( 2. * ( 1. + nu ) );
            kappa  = E / ( 3. * ( 1. - 2.*nu ) );

            // Compute strain Heaviside function
            gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,d,d);

            E_vol = E.trace(); //volumetric strain (scalar) 
            H_pos = (E_vol > 0) ? 1 : 0; // positive Heaviside function
            H_neg = (-E_vol > 0) ? 1 : 0; // negative Heaviside function

            // Compute C (in Voigt)
            result.reshapeCol(i,sz,sz) = math::pow((1. - damage),2) * (kappa*H_pos*Cvol + 2.*mu*Cdev) + kappa*H_neg*Cvol; 
        }

    }

protected:
    using Base::m_data;

};

}
