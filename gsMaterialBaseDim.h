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


#include <gsElasticity/gsMaterialBase.h>

namespace gismo
{

template<class T> 
class gsMaterialData;

template <short_t d, class T>
class gsMaterialBaseDim : public gsMaterialBase<T>
{
public:

    using Base = gsMaterialBase<T>;

    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    /// Shared pointer for gsMaterialBaseDim
    typedef memory::shared_ptr< gsMaterialBaseDim > Ptr;

    /// Unique pointer for gsMaterialBaseDim
    typedef memory::unique_ptr< gsMaterialBaseDim > uPtr;


    gsMaterialBaseDim()
    :
    Base(nullptr,nullptr,nullptr)
    {
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialBaseDim( const gsFunctionSet<T> * mp)
    :
    Base(mp,nullptr,nullptr)
    {
        GISMO_ASSERT(mp->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialBaseDim( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def)
    :
    Base(mp,mp_def,nullptr)
    {
        GISMO_ASSERT(mp->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        GISMO_ASSERT(mp_def->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    gsMaterialBaseDim( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def,
                    const gsFunctionSet<T> * density)
    :
    Base(mp,mp_def,density)
    {
        GISMO_ASSERT(mp->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        GISMO_ASSERT(mp_def->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        this->defaultOptions();
        membersSetZero();
    }

    /// Destructor
    virtual ~gsMaterialBaseDim() {};

public:

    /// See \ref gsMaterialBase for details
    virtual short_t dim() const { return d; }

    /// See \ref gsMaterialBase for details
    virtual void defaultOptions() override
    { }

    /// See \ref gsMaterialBase for details
    virtual inline void eval_deformation_gradient_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const override
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u);

        const index_t N = u.cols();
        gsMatrix<T,d,d> I = gsMatrix<T>::Identity(d,d);
        gsMatrix<T,d,d> physDispJac;
        result.resize(d*d,N);
        for (index_t i=0; i!=N; i++)
        {
            gsAsMatrix<T, Dynamic,Dynamic> jacdef = m_data.mine().m_jac_def.reshapeCol(i,d,d);
            gsAsMatrix<T, Dynamic,Dynamic> jacori = m_data.mine().m_jac_ori.reshapeCol(i,d,d);

            physDispJac = jacdef*(jacori.cramerInverse());
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = result.reshapeCol(i,d,d);
            F = I + physDispJac;
        }
    }

    /// See \ref gsMaterialBase for details
    virtual inline void eval_strain_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const override
    {
        const index_t N = u.cols();
        gsMatrix<T> Fres;
        // NOTE: Could be faster to compute F here in the loop.
        this->eval_deformation_gradient_into(patch,u,Fres);
        gsMatrix<T,d,d> I = gsMatrix<T>::Identity(d,d);
        result.resize(d*d,N);
        for (index_t i=0; i!=N; i++)
        {
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = Fres.reshapeCol(i,d,d);
            // Green-Lagrange strain, E = 0.5*(F'*F-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic,Dynamic> E = result.reshapeCol(i,d,d);
            E = 0.5 * (F.transpose() * F - I);
        }
    }

    virtual void setUndeformed(const gsFunctionSet<T> * undeformed) override
    {
        GISMO_ASSERT(undeformed->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        Base::setUndeformed(f_ptr);
    }

    virtual void setDeformed(const gsFunctionSet<T> * deformed) override
    {
        GISMO_ASSERT(deformed->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        Base::setDeformed(f_ptr);
    }

    virtual void setUndeformed(const function_ptr undeformed) override
    {
        GISMO_ASSERT(undeformed->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        Base::setUndeformed(undeformed);
    }

    virtual void setDeformed(const function_ptr deformed) override
    {
        GISMO_ASSERT(deformed->targetDim()==d,"Geometric dimension and the template dimension are not the same!");
        Base::setDeformed(deformed);
    }

protected:
    void _computeGeometricData(index_t patch, const gsMatrix<T>& u) const
    {
        m_data.mine().m_jac_ori.resize(d*d,u.cols());
        m_data.mine().m_jac_def.resize(d*d,u.cols());

        gsMapData<T> map_ori, map_def;
        map_def.flags = map_ori.flags = NEED_JACOBIAN;
        map_def.points = map_ori.points = u;
        static_cast<const gsFunction<T>&>(Base::m_undeformed->piece(patch)   ).computeMap(map_ori);
        static_cast<const gsFunction<T>&>(Base::m_deformed  ->piece(patch)   ).computeMap(map_def);

        for (index_t k=0; k!= u.cols(); k++)
        {
            m_data.mine().m_jac_ori.reshapeCol(k,d,d) = map_ori.jacobian(k);
            m_data.mine().m_jac_def.reshapeCol(k,d,d) = map_def.jacobian(k);
        }   
    }

    void _computeParameterData(index_t patch, const gsMatrix<T>& u) const
    {
        m_data.mine().m_parmat.resize(m_pars.size(),u.cols());
     
        gsMatrix<T> tmp;

        gsMapData<T> map;
        map.flags = NEED_VALUE;
        map.points = u;
        static_cast<const gsFunction<T>&>(Base::m_undeformed->piece(patch)   ).computeMap(map);

        m_data.mine().m_parmat.resize(m_pars.size(),map.values[0].cols());
        m_data.mine().m_parmat.setZero();

        for (size_t v=0; v!=m_pars.size(); v++)
        {
            m_pars[v]->piece(patch).eval_into(map.values[0], tmp);
            m_data.mine().m_parmat.row(v) = tmp;
        }
    }


protected:

    void membersSetZero()
    {   
        m_data.mine().membersSetZero();
    }

protected:

    using Base::m_undeformed;
    using Base::m_deformed;

    using Base::m_pars;
    using Base::m_density;

    using Base::m_options;

    // Geometric data point
    mutable util::gsThreaded< gsMaterialData<T> > m_data;

public:

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //must be present whenever the class contains fixed size matrices
#   undef Eigen

};


template<class T>
class gsMaterialData
{

public:

    void membersSetZero()
    {
        m_jac_ori.setZero();
        m_jac_def.setZero();
        m_rhoMat.setZero();
    }

    mutable gsMatrix<T> m_parmat;
    mutable gsMatrix<T> m_rhoMat;
    mutable gsMatrix<T> m_jac_ori, m_jac_def;
};

}
