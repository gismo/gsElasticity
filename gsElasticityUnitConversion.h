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

class ElasticityConstantsConverter
{
private:
    real_t _k;
    real_t _e;
    real_t _l;
    real_t _g;
    real_t _n;
    real_t _m;
private:
    static inline real_t abs(const real_t& a) { return  a>=0 ? a : -a; }
public:
    bool operator== (const ElasticityConstantsConverter &other)
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
    const real_t& getK() const
    {
        return _k;
    }
    const real_t& getE() const
    {
        return _e;
    }
    const real_t& getL() const
    {
        return _l;
    }
    const real_t& getG() const
    {
        return _g;
    }
    const real_t& getN() const
    {
        return _n;
    }
    const real_t& getM() const
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
    void setKE(const real_t &k, const real_t &e )
    {
        _k = k;
        _e = e;
        _l = 3*k*(3*k-e)/(9*k-e);
        _g = 3*k*e/(9*k-e);
        _n = (3*k-e)/(6*k);
        _m = 3*k*(3*k+e)/(9*k-e);
    }
    void setKL(const real_t &k, const real_t &l )
    {
        _e = 9*k*(k-l)/(3*k-l);
        setKE(k,_e);
    }
    void setKG(const real_t &k, const real_t &g )
    {
        _e = 9*k*g/(3*k+g);
        setKE(k,_e);
    }
    void setKN(const real_t &k, const real_t &n )
    {
        _e = 3*k*(1-2*n);
        setKE(k,_e);
    }
    void setEG(const real_t &e, const real_t &g )
    {
        _k = e*g/(9*g-3*e);
        setKE(_k,e);
    }
    void setEN(const real_t &e, const real_t &n )
    {
        _k = e/(3-6*n);
        setKE(_k,e);
    }
    void setLG(const real_t &l, const real_t &g )
    {
        _k = l+2*g/3;
        _e = g*(3*l+2*g)/(l+g);
        setKE(_k,_e);

    }
    void setLN(const real_t &l, const real_t &n )
    {
        _k = l*(1+n)/(3*n);
        _e = l*(1+n)*(1-2*n)/n;
        setKE(_k,_e);
    }
    void setGN(const real_t &g, const real_t &n )
    {
        _k = 2*g*(1+n)/(3-6*n);
        _e = 2*g*(1+n);
        setKE(_k,_e);
    }
    void setGM(const real_t &g, const real_t &m )
    {
        _k = m-4*g/3;
        _e = g*(3*m-4*g)/(m-g);
        setKE(_k,_e);
    }
};

std::ostream &operator<< (std::ostream &out, const ElasticityConstantsConverter &data)
{
    out<<"K = "<<data.getK()<<", E = "<<data.getE()<<", L = "<<data.getL()<<", G = "<<data.getG()<<", N = "<<data.getN()<<", M = "<<data.getM()<<"\n";
    return out;
}


}

