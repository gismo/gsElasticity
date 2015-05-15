/** @file gsElasticityUnitConversion.h

    @brief Linear

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#include <gsCore/gsConfig.h>
#pragma once

namespace gismo {



/** @brief Conversion utility for linear elasticity constants
 *  this is not very stable: when tested on 10e6 random values
 *  the relative error can go up till 10e-8
 */
template <typename T=real_t>
class ElasticityConstantsConverter
{
private:
    T _k;
    T _e;
    T _l;
    T _g;
    T _n;
    T _m;
private:
    static inline T abs(const T& a) { return  a>=0 ? a : -a; }
public:
    bool operator== (const ElasticityConstantsConverter<T> &other)
    {
        bool passed = true;
        double tolerance = 10e-6;
        passed = passed && abs(_k-other._k)<abs(_k)*tolerance;
        passed = passed && abs(_e-other._e)<abs(_e)*tolerance;
        passed = passed && abs(_l-other._l)<abs(_l)*tolerance;
        passed = passed && abs(_g-other._g)<abs(_g)*tolerance;
        passed = passed && abs(_n-other._n)<abs(_n)*tolerance;
        passed = passed && abs(_m-other._m)<abs(_m)*tolerance;
        return passed;
    }
    const T& getK() const
    {
        return _k;
    }
    const T& getE() const
    {
        return _e;
    }
    const T& getL() const
    {
        return _l;
    }
    const T& getG() const
    {
        return _g;
    }
    const T& getN() const
    {
        return _n;
    }
    const T& getM() const
    {
        return _m;
    }

    ///
    /// \brief setKE
    /// \param k bulk modulus
    /// \param e Young modulus
    /// compute all the constants from bulk modulus \a k and Young modulus \a e
    /// There are similar function taking different parameters: all possible entries are
    /// \param k bulk modulus
    /// \param e Young modulus
    /// \param l first Lamé coefficient
    /// \param g shear modulus or second Lamé coefficient
    /// \param n Poisson ratio
    /// \param m P-wave modulus
    ///
    ///
    void setKE(const T &k, const T &e )
    {
        _k = k;
        _e = e;
        _l = 3*k*(3*k-e)/(9*k-e);
        _g = 3*k*e/(9*k-e);
        _n = (3*k-e)/(6*k);
        _m = 3*k*(3*k+e)/(9*k-e);
    }
    void setKL(const T &k, const T &l )
    {
        _e = 9*k*(k-l)/(3*k-l);
        setKE(k,_e);
    }
    void setKG(const T &k, const T &g )
    {
        _e = 9*k*g/(3*k+g);
        setKE(k,_e);
    }
    void setKN(const T &k, const T &n )
    {
        _e = 3*k*(1-2*n);
        setKE(k,_e);
    }
    void setEG(const T &e, const T &g )
    {
        _k = e*g/(9*g-3*e);
        setKE(_k,e);
    }
    void setEN(const T &e, const T &n )
    {
        _k = e/(3-6*n);
        setKE(_k,e);
    }
    void setLG(const T &l, const T &g )
    {
        _k = l+2*g/3;
        _e = g*(3*l+2*g)/(l+g);
        setKE(_k,_e);

    }
    void setLN(const T &l, const T &n )
    {
        _k = l*(1+n)/(3*n);
        _e = l*(1+n)*(1-2*n)/n;
        setKE(_k,_e);
    }
    void setGN(const T &g, const T &n )
    {
        _k = 2*g*(1+n)/(3-6*n);
        _e = 2*g*(1+n);
        setKE(_k,_e);
    }
    void setGM(const T &g, const T &m )
    {
        _k = m-4*g/3;
        _e = g*(3*m-4*g)/(m-g);
        setKE(_k,_e);
    }
};

template <typename T>
std::ostream &operator<< (std::ostream &out, const ElasticityConstantsConverter<T> &data)
{
    out<<"K = "<<data.getK()<<", E = "<<data.getE()<<", L = "<<data.getL()<<", G = "<<data.getG()<<", N = "<<data.getN()<<", M = "<<data.getM()<<"\n";
    return out;
}


}

