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

#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsElasticity/gsMaterialUtils.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

template <class T>
class gsLinearDegradedMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsLinearDegradedMaterial(   const T E,
                                const T nu,
                                const gsFunctionSet<T> & damage,
                                short_t d)
    :
    gsLinearDegradedMaterial(gsConstantFunction<T>(E,d),
                             gsConstantFunction<T>(nu,d),damage)
    {
    }

    gsLinearDegradedMaterial(   const gsFunctionSet<T> & E,
                                const gsFunctionSet<T> & nu,
                                const gsFunctionSet<T> & damage)
    :
    Base()
    {
        this->setParameter(0,E);
        this->setParameter(1,nu);
        this->setParameter(2,damage);
    }

    void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Resize the result
        Sresult.resize(dim*dim,N);

        // Lamé parameters
        T E, nu;
        T lambda, G, bulk;
        // Damage variable ("parameter")
        T damage;

        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        gsMatrix<T> E_dev;
        T E_vol, E_vol_pos, E_vol_neg;
        // Loop over points
        for (index_t i=0; i!=N; i++)
        {
            E = data.m_parmat(0,i);
            nu= data.m_parmat(1,i);
            damage = data.m_parmat(2,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            G      = E / ( 2. * ( 1. + nu ) );
            bulk   = lambda + (2./3.)*G;

            gsAsMatrix<T, Dynamic, Dynamic> Emat = data.m_E.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> S = Sresult.reshapeCol(i,dim,dim);

            E_vol = Emat.trace(); //volumetric strain (scalar)
            // E_dev = Emat - 1./dim * E_vol * I; // deviatoric strain
            E_dev = Emat - 1./3 * E_vol * I; // deviatoric strain

            // E_vol_pos = (E_vol >= 0) ? E_vol : 0.; // positive Heaviside function
            // E_vol_neg = (E_vol < 0) ? E_vol : 0.; // negative Heaviside function
            // E_dev = Emat - 1./dim * E_vol * I; // deviatoric strain

            // S = math::pow((1. - damage),2) * (bulk*E_vol_pos*I + 2.*G*E_dev) + bulk*E_vol_neg*I;

            ///// IMPROVE THIS
            T tol = 1e-9;
            T omega = math::pow((1. - damage),2);
            if      (E_vol>tol)
                S = omega*(2.*G*E_dev +  bulk*E_vol*I);
            else if (E_vol<tol)
                S = omega*(2.*G*E_dev) + bulk*E_vol*I;
            else
                S = omega*(2.*G*E_dev);
        }

        // const index_t sz = (dim+1)*dim/2;
        // gsMatrix<T> cmat;
        // eval_matrix_into(data,cmat);
        // gsMatrix<T> Evec;
        // gsMatrix<T> EE= data.m_E.reshapeCol(0,dim,dim);
        // gsMatrix<T> C = cmat.reshapeCol(0,sz,sz);
        // gsMatrix<T> S = Sresult.reshapeCol(0,dim,dim);
        // calculate_voigt_strain(EE,dim,Evec);
        // gsDebugVar((Evec.transpose() * C * Evec).value());
        // gsDebugVar((EE * S ).trace());
        // gsDebugVar(C*Evec);
        // gsDebugVar(S);
        // gsDebugVar(EE);
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
        gsMatrix<T> C_vol = volumetricTensor<T>(dim);
        gsMatrix<T> C_dev = 2.0*deviatoricTensor<T>(dim);
        gsMatrix<T> C_pos;

        // Lamé parameters
        T E, nu;
        T lambda, G, bulk;
        // Damage parameter
        T damage;

        bool H_pos, H_neg;
        T E_vol;

        for (index_t i=0; i!=N; i++)
        {
            E = data.m_parmat(0,i);
            nu= data.m_parmat(1,i);
            damage = data.m_parmat(2,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            G      = E / ( 2. * ( 1. + nu ) );
            bulk   = lambda + (2./3.)*G;

            // Compute strain Heaviside function
            gsAsMatrix<T, Dynamic, Dynamic> Emat = data.m_E.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic, Dynamic> C = Cresult.reshapeCol(i,sz,sz);

            E_vol = Emat.trace(); //volumetric strain (scalar)

            // T tol = math::limits::epsilon();
            // H_pos = (E_vol >= -math::limits::epsilon()) ? 1 : 0; // positive Heaviside function
            // H_neg = (E_vol <  -math::limits::epsilon()) ? 1 : 0; // negative Heaviside function

            // C = math::pow((1. - damage),2) * (bulk*H_pos*C_vol + 2.*G*C_dev) + bulk*H_neg*C_vol;

            ///// IMPROVE THIS
            T tol = 1e-9;
            T omega = math::pow((1. - damage),2);
            if      (E_vol>tol)
                C = omega*(G*C_dev +  bulk*C_vol);
            else if (E_vol<tol)
                C = omega*(G*C_dev) + bulk*C_vol;
            else
                C = omega*(G*C_dev);
        }
    }

    void eval_energy_into(const gsMaterialData<T> & data, gsMatrix<T> & Presult) const
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        // Voigt-size of the tensor
        const index_t sz = (dim+1)*dim/2;

        // Resize the result
        Presult.resize(1,N);

        // Identity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);

        // C tensors
        gsMatrix<T> C_vol = volumetricTensor<T>(dim);
        gsMatrix<T> C_dev = 2.0*deviatoricTensor<T>(dim);
        gsMatrix<T> C_pos;

        // Lamé parameters
        T E, nu;
        T lambda, G, bulk;
        // Damage parameter
        // T damage;

        T E_pos;

        gsMatrix<T> E_vec, tmpE;
        for (index_t i=0; i!=N; i++)
        {
            E = data.m_parmat(0,i);
            nu= data.m_parmat(1,i);
            // damage = data.m_parmat(2,i);
            lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            G      = E / ( 2. * ( 1. + nu ) );
            bulk   = lambda + (2./3.)*G;

            gsAsMatrix<T, Dynamic, Dynamic> Emat = data.m_E.reshapeCol(i,dim,dim);
            tmpE = Emat;
            calculate_voigt_strain(tmpE,dim,E_vec); //E_vect in voigt
            // E_pos = math::max(0.0,Emat.trace()); //volumetric strain (scalar)
            // // H_pos = math::max(0.0,E_pos); // positive Heaviside function
            // // C_pos = bulk*H_pos*C_vol + G*C_dev;
            // C_pos = bulk*E_pos*C_vol + G*C_dev;

            ///// IMPROVE THIS
            T tol = 1e-9;
            if (Emat.trace() > tol)
                C_pos = G*C_dev + bulk*C_vol;
            else
                C_pos = G*C_dev;

            Presult(0,i) = 0.5 * (E_vec.transpose() * C_pos * E_vec).value();
        }
    }
};

}
