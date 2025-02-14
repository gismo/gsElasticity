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

#include <gsCore/gsFunctionSet.h>
#include <gsIO/gsOptionList.h>
#include <gsUtils/gsThreaded.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template<class T>
class gsMaterialData;


template <class T>
class gsMaterialBase
{
public:

    typedef typename gsFunctionSet<T>::Ptr function_ptr;

    /// Shared pointer for gsMaterialBase
    typedef memory::shared_ptr< gsMaterialBase > Ptr;

    /// Unique pointer for gsMaterialBase
    typedef memory::unique_ptr< gsMaterialBase > uPtr;


    gsMaterialBase()
    :
    gsMaterialBase(nullptr,nullptr,nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> * mp)
    :
    gsMaterialBase(mp,nullptr,nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def)
    :
    gsMaterialBase(mp,mp_def,nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> & mp)
    :
    gsMaterialBase(&mp,nullptr,nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> & mp,
                    const gsFunctionSet<T> & mp_def)
    :
    gsMaterialBase(&mp,&mp_def,nullptr)
    {
    }

    // gsMaterialBase( const gsFunctionSet<T> & mp)
    // :
    // gsMaterialBase(mp,nullptr,nullptr)
    // {
    // }

    // gsMaterialBase( const gsFunctionSet<T> & mp,
    //                 const gsFunctionSet<T> & mp_def)
    // :
    // gsMaterialBase(mp,mp_def,nullptr)
    // {
    // }

    gsMaterialBase( const gsFunctionSet<T> & mp,
                    const gsFunctionSet<T> & mp_def,
                    const gsFunctionSet<T> & Density)
    :
    gsMaterialBase(memory::make_shared(mp.clone().release()),
                   memory::make_shared(mp_def.clone().release()),
                   memory::make_shared(Density.clone().release()))
    {
    }

    gsMaterialBase( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def,
                    const gsFunctionSet<T> * Density)
    :
    gsMaterialBase(memory::make_shared_not_owned(mp),
                   memory::make_shared_not_owned(mp_def),
                   memory::make_shared_not_owned(Density))
    {
    }

    gsMaterialBase( const function_ptr & mp,
                    const function_ptr & mp_def,
                    const function_ptr & Density)
    :
    m_density(Density)
    {
        m_undeformed.mine() = mp;
        m_deformed.mine() = mp_def;
    }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsMaterialBase, clone)

    /// Destructor
    virtual ~gsMaterialBase() {};

    /// Copy constructor (makes deep copy)
    gsMaterialBase( const gsMaterialBase<T> & other )
    {
        operator=(other);
    }

    /// Move constructor
    gsMaterialBase( gsMaterialBase<T> && other )
    {
        operator=(give(other));
    }

    gsMaterialBase<T> & operator= ( const gsMaterialBase<T> & other )
    {
        if (this!=&other)
        {
            m_undeformed.mine() = other.m_undeformed.mine();
            m_deformed.mine() = other.m_deformed.mine();
            m_options = other.m_options;

            m_pars = other.m_pars;
            m_density = other.m_density;
        }
        return *this;
    }

    gsMaterialBase<T> & operator= ( gsMaterialBase<T> && other )
    {
        m_undeformed.mine() = give(other.m_undeformed.mine());
        m_deformed.mine() = give(other.m_deformed.mine());
        m_options = give(other.m_options);


        m_pars = give(other.m_pars);
        m_density = give(other.m_density);
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
    virtual void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /**
     *
     */
    virtual void precompute(const index_t patch, const gsMatrix<T> & u) const
    {
        this->_computeGeometricData(patch,u,m_data.mine());
    }

    /**
     * @brief      Evaluates the deformation gradient at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_deformation_gradient_into(gsMatrix<T> & Fresult)
    {
        this->eval_deformation_gradient_into(m_data.mine(),Fresult);
    }

    virtual void eval_deformation_gradient_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Fresult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_deformation_gradient_into(m_data.mine(),Fresult);
    }

    virtual void eval_deformation_gradient_into(const gsMaterialData<T> & data, gsMatrix<T> & Fresult) const
    {
        Fresult = data.m_F;
    }


    /**
     * @brief      Evaluates the strain tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_strain_into(gsMatrix<T> & Eresult)
    {
        this->eval_strain_into(m_data.mine(),Eresult);
    }

    virtual void eval_strain_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Eresult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_strain_into(m_data.mine(),Eresult);
    }

    virtual void eval_strain_into(const gsMaterialData<T> & data, gsMatrix<T> & Eresult) const
    {
        Eresult = data.m_E;
    }


    /**
     * @brief      Evaluates the stress tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_stress_into(gsMatrix<T> & Sresult)
    {
        this->eval_stress_into(m_data.mine(),Sresult);
    }

    virtual void eval_stress_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Sresult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_stress_into(m_data.mine(),Sresult);
    }

    virtual void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the material tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_matrix_into(gsMatrix<T> & Cresult)
    {
        this->eval_matrix_into(m_data.mine(),Cresult);
    }

    virtual void eval_matrix_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Cresult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_matrix_into(m_data.mine(),Cresult);
    }

    virtual void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the undamaged positive material tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_matrix_pos_into(gsMatrix<T> & Cresult)
    {
        this->eval_matrix_pos_into(m_data.mine(),Cresult);
    }

    virtual void eval_matrix_pos_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Cresult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_matrix_pos_into(m_data.mine(),Cresult);
    }

    virtual void eval_matrix_pos_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the energy at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_energy_into(gsMatrix<T> & Presult)
    {
        this->eval_energy_into(m_data.mine(),Presult);
    }

    virtual void eval_energy_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & Presult) const
    {
        // Compute map and parameters
        this->_computeGeometricData(patch,u,m_data.mine());
        this->eval_energy_into(m_data.mine(),Presult);
    }

    virtual void eval_energy_into(const gsMaterialData<T> & data, gsMatrix<T> & Presult) const
    { GISMO_NO_IMPLEMENTATION; }

    /// Sets the density
    virtual inline void setDensity(function_ptr Density) { m_density = Density; }
    /// Sets the density
    virtual inline void setDensity(const gsFunctionSet<T> & Density)
    {
        function_ptr fun = memory::make_shared(Density.clone().release());
        m_density = fun;
    }
    /// Returns true if a density is assigned
    virtual inline bool hasDensity() const { return m_density!=nullptr; }

    /// Gets the Density
    virtual const function_ptr getDensity()  const {return m_density;}

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual inline void setParameters(const std::vector<function_ptr> &pars)
    {
        m_pars = pars;
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual inline void setParameter(const index_t i, const function_ptr &par)
    {
        m_pars[i] = par;
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual inline void setParameters(const std::vector<gsFunctionSet<T> *> &pars)
    {
        m_pars.resize(pars.size());
        for (size_t k = 0; k!=pars.size(); k++)
            m_pars[k] = memory::make_shared_not_owned(pars[k]);
    }

    /**
     * @brief      Sets the material parameters.
     *
     * @param[in]  pars  Function pointers for the parameters in a container
     */
    virtual inline void setParameter(const index_t i, const gsFunctionSet<T> &par)
    {
        if ((index_t)m_pars.size() < i+1)
            m_pars.resize(i+1);
        m_pars[i] = memory::make_shared(par.clone().release());
    }

    /**
     * @brief      Gets parameter i
     *
     * @param[in]  i     The parameter index
     *
     * @return     The parameter.
     */
    virtual inline const function_ptr getParameter(const index_t i)  const
    {
        GISMO_ASSERT(i < (index_t)m_pars.size(),"Parameter "<<i<<" is unavailable");
        return m_pars[i] ;
    }

    /**
     * @brief      Gets the number of parameters
     *
     */
    virtual inline index_t numParameters() const { return m_pars.size(); }

    /// See \ref gsMaterialMatrixBase for details
    virtual inline void resetParameters()
    {
        m_pars.clear();
        m_pars.resize(0);
    }

    virtual inline void setUndeformed(const gsFunctionSet<T> * undeformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        m_undeformed.mine() = f_ptr;
    }
    virtual inline void setDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        m_deformed.mine() = f_ptr;
    }

    virtual inline void setUndeformed(const function_ptr undeformed)
    {
        m_undeformed.mine() = undeformed;
    }
    virtual inline void setDeformed(const function_ptr deformed)
    {
        m_deformed.mine() = deformed;
    }

    const inline function_ptr getUndeformed() const { return m_undeformed.mine(); }
    const inline function_ptr getDeformed()  const { return m_deformed.mine(); }

    virtual inline bool hasUndeformed() const
    { return m_deformed.mine()!=nullptr; }
    virtual inline bool hasDeformed() const { return m_deformed.mine()!=nullptr; }

    void precomputeData(index_t patch, const gsMatrix<T>& u)
    {
        this->computeGeometricData(patch,u);
    }

    /**
     *
     */
    const gsMaterialData<T> & getData() const{ return m_data.mine(); }

    void setData(const gsMaterialData<T> & data) { m_data.mine() = data; }

protected:

    GISMO_DEPRECATED
    void _computeGeometricData(index_t patch, const gsMatrix<T>& u) const
    {
        this->_computeGeometricData(patch,u,m_data.mine());
    }

    void _computeGeometricData(index_t patch, const gsMatrix<T>& u, gsMaterialData<T> & data) const
    {
        data.dim = u.rows();
        data.size = u.cols();
        GISMO_ASSERT(m_undeformed.mine()!=nullptr,"Undeformed function is not set");
        GISMO_ASSERT(m_deformed.mine()!=nullptr,"Deformed function is not set");
        GISMO_ASSERT(m_undeformed.mine()->domainDim()==data.dim,"Geometric dimension and the point dimension are not the same!");
        GISMO_ASSERT(m_deformed.mine()->domainDim()==data.dim,"Geometric dimension and the point dimension are not the same!");
        data.m_F.resize(data.dim*data.dim,data.size);
        data.m_E.resize(data.dim*data.dim,data.size);
        data.m_parmat.resize(m_pars.size(),u.cols());

        gsMapData<T> map_ori, map_def;
        map_ori.flags = NEED_VALUE;
        map_ori.flags|= NEED_JACOBIAN;
        map_def.flags = NEED_JACOBIAN;
        map_def.points = map_ori.points = u;
        m_undeformed.mine()->function(patch).computeMap(map_ori);
        m_deformed.mine()->function(patch).computeMap(map_def);

        data.m_parmat.setZero();
        gsMatrix<T> pars;
        gsMatrix<T> jac_ori, jac_def;
        gsMatrix<T> I = gsMatrix<T>::Identity(data.dim,data.dim);
        for (index_t k=0; k!= data.size; k++)
        {
            jac_ori = map_ori.jacobian(k);
            jac_def = map_def.jacobian(k);
            // deformation gradient F = I + du/dx
            gsAsMatrix<T, Dynamic,Dynamic> F = data.m_F.reshapeCol(k,data.dim,data.dim);
            F = jac_def*(jac_ori.cramerInverse());;
            // strain tensor E = 0.5*(F'*F-I), a.k.a. full geometric strain tensor
            gsAsMatrix<T, Dynamic,Dynamic> E = data.m_E.reshapeCol(k,data.dim,data.dim);
            E = 0.5 * (F.transpose() * F - I);
        }

        for (size_t v=0; v!=m_pars.size(); v++)
        {
            m_pars[v]->piece(patch).eval_into(map_ori.values[0], pars);
            data.m_parmat.row(v) = pars;
        }
    }

    // void _computeParameterData(index_t patch, const gsMatrix<T>& u) const
    // {
    //     this->_computeParameterData(patch,u,m_data.mine());
    // }

    // void _computeParameterData(index_t patch, const gsMatrix<T>& u, gsMaterialData<T> & data) const
    // {

    //     gsMatrix<T> tmp;

    //     gsMapData<T> map;
    //     map.flags = NEED_VALUE;
    //     map.points = u;
    //     static_cast<const gsFunction<T>&>(m_undeformed.mine()->piece(patch)   ).computeMap(map);

    //     data.m_parmat.resize(m_pars.size(),map.values[0].cols());
    //     data.m_parmat.setZero();

    //     for (size_t v=0; v!=m_pars.size(); v++)
    //     {
    //         m_pars[v]->piece(patch).eval_into(map.values[0], tmp);
    //         data.m_parmat.row(v) = tmp;
    //     }
    // }


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
        m_F.setZero();
        m_E.setZero();
        m_rhoMat.setZero();
    }

    mutable gsMatrix<T> m_parmat;
    mutable gsMatrix<T> m_rhoMat;
    // mutable gsMatrix<T> m_jac_ori, m_jac_def;
    mutable gsMatrix<T> m_F, m_E;

    mutable short_t dim;
    mutable index_t size;
};


}
