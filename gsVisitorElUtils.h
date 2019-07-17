/** @file gsVisitorElUtils.h

    @brief Several special matrix forming functions for elasticity visitors.

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

// construct an elasticity tensor C in Voigt notation assuming that
// C = a R x R + b R (*) R in tensor notation,
// where R is a symmetric second order tensoR (letter T is reserved).
// (see Bernal, Calo, Collier, et. at., "PetIGA", ICCS 2013, p. 1608)
template <class T>
inline void setC(gsMatrix<T> & C, const gsMatrix<T> & R, T a, T b)
{
    short_t dim = R.cols();
    short_t dimTensor = dim*(dim+1)/2;
    C.resize(dimTensor,dimTensor);
    // voigt indices
    gsMatrix<unsigned> v(dimTensor,2);
    if (dim == 2)
        v << 0,0,
             1,1,
             0,1;
    if (dim == 3)
        v << 0,0,
             1,1,
             2,2,
             0,1,
             1,2,
             0,2;
    // componentwise definition
    for (short_t i = 0; i < dimTensor; ++i)
        for (short_t j = 0; j < dimTensor; ++j)
            C(i,j) = a*R(v(i,0),v(i,1))*R(v(j,0),v(j,1)) +
                     b*(R(v(i,0),v(j,0))*R(v(i,1),v(j,1)) +
                        R(v(i,0),v(j,1))*R(v(i,1),v(j,0)));
}

// transform stress tensor S to a vector in Voigt notation
template <class T>
inline void voigtStress(gsVector<T> & Svec, const gsMatrix<T> & S)
{
    short_t dim = S.cols();
    short_t dimTensor = dim*(dim+1)/2;
    Svec.resize(dimTensor);
    if (dim == 2)
    {
        Svec(0) = S(0,0);
        Svec(1) = S(1,1);
        Svec(2) = S(0,1);
    }
    if (dim == 3)
    {
        Svec(0) = S(0,0);
        Svec(1) = S(1,1);
        Svec(2) = S(2,2);
        Svec(3) = S(0,1);
        Svec(4) = S(1,2);
        Svec(5) = S(0,2);
    }
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
