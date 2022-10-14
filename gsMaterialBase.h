/** @file gsMaterialBase.h

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

// #include <gsUtils/gsThreaded.h>

namespace gismo
{

template <class T>
class gsMaterialBase
{

public:
    virtual ~gsMaterialBase() {};

    gsMaterialBase( const gsFunctionSet<T> * undeformed,
                    const gsFunctionSet<T> * deformed)
    :
    m_undeformed(undeformed),
    m_deformed(deformed)
    {

    }

    gsMaterialBase( const gsFunctionSet<T> * undeformed)
    :
    m_undeformed(undeformed),
    m_deformed(nullptr)
    {

    }

    /**
     * @brief      { function_description }
     *
     * @param[in]  u       { parameter_description }
     * @param      result  The result
     */
    virtual void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result,
                                  index_t k = 0) const {GISMO_NO_IMPLEMENTATION};// = 0;

    /**
     * @brief      { function_description }
     *
     * @param[in]  u       { parameter_description }
     * @param      result  The result
     */
    virtual void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result,
                                  index_t k = 0) const {GISMO_NO_IMPLEMENTATION};// = 0;

    void gradient_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Fres, index_t k = 0) const
    {
        gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
        gsMapData<T> map_def = _getMapData(deformed,u,k);

        gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
        gsMatrix<T> physDispJac;
        Fres.resize(u.rows()*u.rows(),u.cols());
        for (index_t i=0; i!=u.cols(); i++)
        {
            physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,u.rows(),u.rows());
            F = I + physDispJac;
        }
    }

    void gradient_into(const gsMatrix<T>& u, gsMatrix<T>& Fres, index_t k = 0) const
    {
        this->gradient_into(m_deformed,u,Fres,k);
    }

    void strain_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Eres, index_t k = 0) const
    {
        gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
        gsMapData<T> map_def = _getMapData(deformed,u,k);

        gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
        gsMatrix<T> physDispJac, F, RCG;
        Eres.resize(u.rows()*u.rows(),u.cols());
        for (index_t i=0; i!=u.cols(); i++)
        {
            physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());

            // deformation gradient F = I + du/dx = dx/dX
            F = I + physDispJac;
            // Right Cauchy Green strain, C = F'*F
            RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
            E = 0.5 * (RCG - I);
        }
    }

    void strain_into(const gsMatrix<T>& u, gsMatrix<T>& Eres, index_t k = 0) const
    {
        this->strain_into(m_deformed,u,Eres,k);
    }

    void allStrains_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Fres, gsMatrix<T>& Eres, index_t k = 0) const
    {
        gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
        gsMapData<T> map_def = _getMapData(deformed,u,k);

        gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
        gsMatrix<T> physDispJac, F, RCG, E;
        Fres.resize(u.rows()*u.rows(),u.cols());
        Eres.resize(u.rows()*u.rows(),u.cols());
        for (index_t i=0; i!=u.cols(); i++)
        {
            physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());
            // deformation gradient F = I + du/dx = dx/dX
            gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,u.rows(),u.rows());
            F = I + physDispJac;
            // Right Cauchy Green strain, C = F'*F
            RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
            E = 0.5 * (RCG - I);
        }
    }

    void allStrains_into(const gsMatrix<T>& u, gsMatrix<T>& Fres, gsMatrix<T>& Eres, index_t k = 0) const
    {
        this->allStrains_into(m_deformed,u,Fres,Eres,k);
    }

    void setDeformed(const gsFunctionSet<T> * deformed)
    {
        m_deformed = deformed;
    }

protected:
    gsMapData<T> _getMapData(const gsFunctionSet<T> * geometry, const gsMatrix<T>& u, index_t k = 0) const
    {
        gsMapData<T> res;
        res.flags  = NEED_JACOBIAN;
        res.points = u;
        static_cast<const gsFunction<T>&>(geometry ->piece(k)).computeMap(res);
        return res;
    }

protected:
    const gsFunctionSet<T> * m_undeformed;
    const gsFunctionSet<T> * m_deformed;


    // mutable util::gsThreaded<gsMapData<T> > map_ori, map_def;

};

}
