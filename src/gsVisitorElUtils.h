/** @file gsVisitorElUtils.h

    @brief Tensor operations for elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

namespace gismo
{

// Indices of the Voigt notation
inline short_t voigt(short_t dim, short_t I, short_t J)
{
    if (dim == 2)
        switch(I)
        {
        case 0: return J == 0 ? 0 : 0;
        case 1: return J == 0 ? 1 : 1;
        case 2: return J == 0 ? 0 : 1;
        }
    else if (dim == 3)
        switch (I)
        {
        case 0: return J == 0 ? 0 : 0;
        case 1: return J == 0 ? 1 : 1;
        case 2: return J == 0 ? 2 : 2;
        case 3: return J == 0 ? 0 : 1;
        case 4: return J == 0 ? 1 : 2;
        case 5: return J == 0 ? 0 : 2;
        }
    GISMO_ERROR("voigt notation indices error");
}

// construct a fourth order symmetric identity tensor C based on a second order symmetric tensor R
// C_ijkl = (R_ik*R_jl + R_il*R_jk)/2 in Voigt notation
template <class T>
inline void symmetricIdentityTensor(gsMatrix<T> & C, const gsMatrix<T> & R)
{
    short_t dim = R.cols();
    short_t dimTensor = dim*(dim+1)/2;
    C.setZero(dimTensor,dimTensor);
    // componentwise definition
    for (short_t i = 0; i < dimTensor; ++i)
        for (short_t j = 0; j < dimTensor; ++j)
            C(i,j) = (R(voigt(dim,i,0),voigt(dim,j,0))*R(voigt(dim,i,1),voigt(dim,j,1)) +
                        R(voigt(dim,i,0),voigt(dim,j,1))*R(voigt(dim,i,1),voigt(dim,j,0)));
}

// Deviatoric matrix operator in 2D and 3D
template <class T> //template it with dimension
inline void deviatoricTensor(gsMatrix<T> & C, const gsMatrix<T> & R)
{
    short_t dim = R.cols();
    short_t dimTensor = dim*(dim+1)/2;
    C.setZero(dimTensor,dimTensor);

    C.block(0,0,dim,dim).setConstant(-1.0/3.0);
    C.block(0,0,dim,dim).diagonal().setConstant(2.0/3.0);
    C.block(dim,dim,dimTensor-dim,dimTensor-dim).diagonal().setConstant(0.5);

    // if (dim == 2) 
    // {
    // C <<2.0/3.0, -1.0/3.0, 0.,
    //     -1.0/3.0, 2.0/3.0, 0.,
    //     0., 0., 1.0/2.0;
    // } 
    // else if (dim == 3) 
    // {
    // C << 2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
    //     -1.0/3.0, 2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
    //     -1.0/3.0, -1.0/3.0, 2.0/3.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 1.0/2.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 1.0/2.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 0.0, 1.0/2.0;
    // }
}


// construct a fourth order matrix-trace tensor C based on two second order symmetric tensors R and S
// C_ijkl = R_ij*S_kl in Voigt notation
template <class T>
inline void matrixTraceTensor(gsMatrix<T> & C, const gsMatrix<T> & R, const gsMatrix<T> & S)
{
    short_t dim = R.cols();
    short_t dimTensor = dim*(dim+1)/2;
    C.setZero(dimTensor,dimTensor);
    // componentwise definition
    for (short_t i = 0; i < dimTensor; ++i)
        for (short_t j = 0; j < dimTensor; ++j)
            C(i,j) = R(voigt(dim,i,0),voigt(dim,i,1))*S(voigt(dim,j,0),voigt(dim,j,1));
}

// transform stress tensor S to a vector in Voigt notation
template <class T>
inline void voigtStress(gsVector<T> & Svec, const gsMatrix<T> & S)
{
    short_t dim = S.cols();
    short_t dimTensor = dim*(dim+1)/2;
    Svec.resize(dimTensor);
    for (short i = 0; i < dimTensor; ++i)
        Svec(i) = S(voigt(dim,i,0),voigt(dim,i,1));
}

// transform stress Svec in Voigt notation to stress tensor S 
template <class T>
inline void tensorStress(const index_t dim, const gsMatrix<T> & Svec, gsMatrix<T> & S)
{
    if (dim == 2)
    {
        S.setZero(dim,dim);
        S(0,0) = Svec(0,0);
        S(1,1) = Svec(1,0);
        S(0,1) = Svec(2,0);
        S(1,0) = Svec(2,0);
    }
    else if (dim == 3)
    {
        S.setZero(dim,dim);
        S(0,0) = Svec(0,0);
        S(1,1) = Svec(1,0);
        S(2,2) = Svec(2,0);

        S(0,1) = Svec(3,0); // tao_xy 
        S(1,0) = Svec(3,0); // tao_xy 

        S(0,2) = Svec(4,0); // tao_yz 
        S(2,0) = Svec(4,0); // tao_yz

        S(1,2) = Svec(5,0); // tao_xz
        S(2,1) = Svec(5,0); // tao_xz
    }
}

// transform strain tensor E to a vector in Voigt notation
template <class T>
inline void voigtStrain(gsVector<T> & Evec, const gsMatrix<T> & E)
{
    short_t dim = E.cols();
    //gsDebugVar(dim);
    short_t dimTensor = dim*(dim+1)/2;
    Evec.resize(dimTensor);
    for (short i = 0; i < dimTensor; ++i)
        if (voigt(dim,i,0) != voigt(dim,i,1))
            Evec(i) = E(voigt(dim,i,0),voigt(dim,i,1)) + E(voigt(dim,i,1),voigt(dim,i,0)); // off-diagonal terms
        else
            Evec(i) = E(voigt(dim,i,0),voigt(dim,i,1)); // diagonal terms
}

// auxiliary matrix B such that E:S = B*Svec in the weak form
// (see Bernal, Calo, Collier, et. at., "PetIGA", ICCS 2013, p. 1610)
template <class T>
inline void setB( gsMatrix<T> & B, const gsMatrix<T> & F, const gsVector<T> & bGrad)
{
    short_t dim = F.cols();
    short_t dimTensor = dim*(dim+1)/2;
    B.resize(dimTensor,dim);

    for (short_t j = 0; j < dim; ++j)
    {
        for (short_t i = 0; i < dim; ++i)
            B(i,j) = F(j,i) * bGrad(i);
        if (dim == 2)
            B(2,j) = F(j,0) * bGrad(1) + F(j,1) * bGrad(0);
        if (dim == 3)
            for (short_t i = 0; i < dim; ++i)
            {
                short_t k = (i+1)%dim;
                B(i+dim,j) = F(j,i) * bGrad(k) + F(j,k) * bGrad(i);
            }
    }
}

} // namespace gismo
