/** @file gsMaterialHelper.h

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

#include <gsCore/gsFunctionSet.h>
#include <gsIO/gsOptionList.h>
#include <gsUtils/gsThreaded.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template<class T>
class gsMaterialData;


template <class T>
class gsMaterialHelper
{
public:

    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    /// Shared pointer for gsMaterialHelper
    typedef memory::shared_ptr< gsMaterialHelper > Ptr;

    /// Unique pointer for gsMaterialHelper
    typedef memory::unique_ptr< gsMaterialHelper > uPtr;


    gsMaterialHelper()
    :
    gsMaterialHelper(nullptr,nullptr)
    {
    }

    gsMaterialHelper( const gsFunctionSet<T> * mp)
    :
    gsMaterialHelper(mp,nullptr)
    {
    }

    gsMaterialHelper( const gsFunctionSet<T> & mp)
    :
    gsMaterialHelper(&mp,nullptr)
    {
    }

    gsMaterialHelper( const gsFunctionSet<T> & mp,
                    const gsFunctionSet<T> & mp_def)
    :
    gsMaterialHelper(memory::make_shared(mp.clone().release()),
                   memory::make_shared(mp_def.clone().release()))
    {
    }

    gsMaterialHelper( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def)
    :
    gsMaterialHelper(memory::make_shared_not_owned(mp),
                   memory::make_shared_not_owned(mp_def))
    {
    }

    gsMaterialHelper( const function_ptr & mp,
                    const function_ptr & mp_def)
    {
        m_undeformed.mine() = mp;
        m_deformed.mine() = mp_def;
    }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsMaterialHelper, clone)

    /// Destructor
    virtual ~gsMaterialHelper() {};

    /// Copy constructor (makes deep copy)
    gsMaterialHelper( const gsMaterialHelper<T> & other )
    {
        operator=(other);
    }

    /// Move constructor
    gsMaterialHelper( gsMaterialHelper<T> && other )
    {
        operator=(give(other));
    }

    gsMaterialHelper<T> & operator= ( const gsMaterialHelper<T> & other )
    {
        if (this!=&other)
        {
            m_undeformed.mine() = other.m_undeformed.mine();
            m_deformed.mine() = other.m_deformed.mine();
            m_options = other.m_options;
        }
        return *this;
    }

    gsMaterialHelper<T> & operator= ( gsMaterialHelper<T> && other )
    {
        m_undeformed.mine() = give(other.m_undeformed.mine());
        m_deformed.mine() = give(other.m_deformed.mine());
        m_options = give(other.m_options);
        return *this;
    }

    /**
     * @brief      Sets the default options
     */
    virtual inline void defaultOptions()
    { }

    /**
     * @brief      Returns the options
     *
     * @return     \ref gsOptionList with the class options
     */
    virtual inline gsOptionList & options() { return m_options; }

    /**
     * @brief      Returns the dimension
     *
     * @return     The dimension
     */
    virtual inline short_t dim() const
    {
        GISMO_ASSERT(m_undeformed.mine()!=nullptr,"Undeformed geometry is not set!");
        return m_undeformed.mine()->domainDim();
    }

    /**
     * @brief      Sets the options
     *
     * @param[in]  opt   \ref gsOptionList
     */
    virtual inline void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /**
     * @brief      Evaluates the deformation gradient at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual inline void eval_deformation_gradient_and_strain_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Fresult, gsMatrix<T> & Eresult) const
    {
        const short_t dim = u.rows();
        // Compute map and parameters
        this->_computeGeometricData(patch,u);

        const index_t N = u.cols();
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        gsMatrix<T> physDefJac;
        Fresult.resize(dim*dim,N);
        Eresult.resize(dim*dim,N);
        for (index_t i=0; i!=N; i++)
        {
            gsAsMatrix<T, Dynamic,Dynamic> jacdef = m_data.mine().m_jac_def.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic,Dynamic> jacori = m_data.mine().m_jac_ori.reshapeCol(i,dim,dim);

            physDefJac = jacdef*(jacori.cramerInverse());
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = Fresult.reshapeCol(i,dim,dim);
            F = physDefJac;

            // strain tensor E = 0.5*(F'*F-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic,Dynamic> E = Eresult.reshapeCol(i,dim,dim);
            E = 0.5 * (F.transpose() * F - I);
        }
    }

    /**
     * @brief      Evaluates the deformation gradient at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual inline void eval_deformation_gradient_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    {
        const short_t dim = u.rows();
        // Compute map and parameters
        this->_computeGeometricData(patch,u);

        const index_t N = u.cols();
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        gsMatrix<T> physDefJac;
        result.resize(dim*dim,N);
        for (index_t i=0; i!=N; i++)
        {
            gsAsMatrix<T, Dynamic,Dynamic> jacdef = m_data.mine().m_jac_def.reshapeCol(i,dim,dim);
            gsAsMatrix<T, Dynamic,Dynamic> jacori = m_data.mine().m_jac_ori.reshapeCol(i,dim,dim);

            physDefJac = jacdef*(jacori.cramerInverse());
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = result.reshapeCol(i,dim,dim);

            F = physDefJac;

        }
    }


    /**
     * @brief      Evaluates the strain tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual inline void eval_strain_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    {
        const short_t dim = u.rows();
        const index_t N = u.cols();
        gsMatrix<T> Fres;
        // NOTE: Could be faster to compute F here in the loop.
        this->eval_deformation_gradient_into(patch,u,Fres);
        gsMatrix<T> I = gsMatrix<T>::Identity(dim,dim);
        result.resize(dim*dim,N);

        for (index_t i=0; i!=N; i++)
        {
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = Fres.reshapeCol(i,dim,dim);
            // Green-Lagrange strain, E = 0.5*(F'*F-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic,Dynamic> E = result.reshapeCol(i,dim,dim);

            E = 0.5 * (F.transpose() * F - I);
            // E = 0.5 * ((F-I).transpose() + (F-I)); // small strains
        }
    }

    /// Returns true if a density is assigned
    virtual inline bool hasDensity() const { return m_density!=nullptr; }

    virtual void setUndeformed(const gsFunctionSet<T> * undeformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        m_undeformed.mine() = f_ptr;
    }
    virtual void setDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        m_deformed.mine() = f_ptr;
    }

    virtual void setUndeformed(const function_ptr undeformed)
    {
        m_undeformed.mine() = undeformed;
    }
    virtual void setDeformed(const function_ptr deformed)
    {
        m_deformed.mine() = deformed;
    }

    const function_ptr getUndeformed() const { return m_undeformed.mine(); }
    const function_ptr getDeformed()  const { return m_deformed.mine(); }

    virtual bool hasUndeformed() const
    { return m_deformed.mine()!=nullptr; }
    virtual bool hasDeformed() const { return m_deformed.mine()!=nullptr; }

    void precompute(index_t patch, const gsMatrix<T>& u)
    {
        this->computeGeometricData(patch,u);
    }

    const gsMaterialData<T> & data() const{ return m_data.mine(); }


    /**
     *
     */
    const gsMaterialData<T> & getData() const{ return m_data.mine(); }

    void setData(const gsMaterialData<T> & data) { m_data.mine() = data; }

protected:
    void _computeGeometricData(index_t patch, const gsMatrix<T>& u) const
    {
        const short_t dim = u.rows();
        GISMO_ASSERT(m_undeformed.mine()!=nullptr,"Undeformed function is not set");
        GISMO_ASSERT(m_deformed.mine()!=nullptr,"Deformed function is not set");
        GISMO_ASSERT(m_undeformed.mine()->domainDim()==dim,"Geometric dimension and the point dimension are not the same!");
        GISMO_ASSERT(m_deformed.mine()->domainDim()==dim,"Geometric dimension and the point dimension are not the same!");
        m_data.mine().m_jac_ori.resize(dim*dim,u.cols());
        m_data.mine().m_jac_def.resize(dim*dim,u.cols());

        gsMapData<T> map_ori, map_def;
        map_def.flags = map_ori.flags = NEED_JACOBIAN;
        map_def.points = map_ori.points = u;
        m_undeformed.mine()->function(patch).computeMap(map_ori);
        m_deformed.mine()->function(patch).computeMap(map_def);

        for (index_t k=0; k!= u.cols(); k++)
        {
            m_data.mine().m_jac_ori.reshapeCol(k,dim,dim) = map_ori.jacobian(k);
            m_data.mine().m_jac_def.reshapeCol(k,dim,dim) = map_def.jacobian(k);
        }
    }

    void _computeParameterData(index_t patch, const gsMatrix<T>& u) const
    {
        const short_t dim = u.rows();
        GISMO_ASSERT(m_undeformed.mine()!=nullptr,"Undeformed function is not set");
        GISMO_ASSERT(m_undeformed.mine()->domainDim()==dim,"Geometric dimension and the point dimension are not the same!");
        m_data.mine().m_parmat.resize(m_pars.size(),u.cols());

        gsMatrix<T> tmp;

        gsMapData<T> map;
        map.flags = NEED_VALUE;
        map.points = u;
        static_cast<const gsFunction<T>&>(m_undeformed.mine()->piece(patch)   ).computeMap(map);

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

    util::gsThreaded<function_ptr> m_undeformed;
    util::gsThreaded<function_ptr> m_deformed;

    std::vector< function_ptr > m_pars;
    function_ptr m_density;

    gsOptionList m_options;

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
