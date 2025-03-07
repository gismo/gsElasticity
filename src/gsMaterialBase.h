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

/*
    NOTE ON THREADING:
    The gsMaterialEvalSingle and gsMaterialEval classes contain a thread-safe data object m_data
    which is used to store the material data for each thread. This is done to avoid

 */

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
    gsMaterialBase(nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> & Density)
    :
    gsMaterialBase(memory::make_shared(Density.clone().release()))
    {
    }

    gsMaterialBase( const gsFunctionSet<T> * Density)
    :
    gsMaterialBase(memory::make_shared_not_owned(Density))
    {
    }

    gsMaterialBase( const function_ptr & Density)
    :
    m_density(Density)
    {
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
            m_options = other.m_options;

            m_pars = other.m_pars;
            m_density = other.m_density;
        }
        return *this;
    }

    gsMaterialBase<T> & operator= ( gsMaterialBase<T> && other )
    {
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
     * @brief      Sets the options
     *
     * @param[in]  opt   \ref gsOptionList
     */
    virtual void setOptions(gsOptionList opt) {m_options.update(opt,gsOptionList::addIfUnknown); }

    /**
     *
     */

    virtual void precompute(const gsFunctionSet<T> * undeformed,
                            const index_t patch,
                            const gsMatrix<T> & u,
                            gsMaterialData<T> & data) const
    {
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(undeformed,nullptr,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
    }

    virtual void precompute(const gsFunctionSet<T> * undeformed,
                            const gsFunctionSet<T> * deformed,
                            const index_t patch,
                            const gsMatrix<T> & u,
                            gsMaterialData<T> & data) const
    {
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(undeformed,deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
    }

    virtual void precompute(const gsMapData<T> & map_ori,
                            const gsMapData<T> & map_def,
                            gsMaterialData<T> & data) const
    {
        this->_computeStrains(map_ori,map_def,data);
    }

    /**
     * @brief      Evaluates the deformation gradient at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return
     */
    virtual void eval_deformation_gradient_into(const gsFunctionSet<T> & undeformed,
                                                const gsFunctionSet<T> & deformed,
                                                const index_t patch,
                                                const gsMatrix<T>& u,
                                                gsMatrix<T> & Fresult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_deformation_gradient_into(data,Fresult);
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
    virtual void eval_strain_into(  const gsFunctionSet<T> & undeformed,
                                    const gsFunctionSet<T> & deformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                    gsMatrix<T> & Eresult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_strain_into(data,Eresult);
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
    virtual void eval_stress_into(  const gsFunctionSet<T> & undeformed,
                                    const gsFunctionSet<T> & deformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                    gsMatrix<T> & Sresult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_stress_into(data,Sresult);
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
    virtual void eval_matrix_into(  const gsFunctionSet<T> & undeformed,
                                    const gsFunctionSet<T> & deformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                    gsMatrix<T> & Cresult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_matrix_into(data,Cresult);
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
    virtual void eval_matrix_pos_into(  const gsFunctionSet<T> & undeformed,
                                        const gsFunctionSet<T> & deformed,
                                        const index_t patch,
                                        const gsMatrix<T>& u,
                                        gsMatrix<T> & Cresult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_matrix_pos_into(data,Cresult);
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
    virtual void eval_energy_into(  const gsFunctionSet<T> & undeformed,
                                    const gsFunctionSet<T> & deformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                    gsMatrix<T> & Presult) const
    {
        // Compute map and parameters
        gsMaterialData<T> data;
        gsMapData<T> map_ori, map_def;
        this->_computeGeometricData(&undeformed,&deformed,patch,u,map_ori,map_def);
        this->_computeStrains(map_ori,map_def,data);
        this->eval_energy_into(data,Presult);
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

protected:

    void _computeStrains(const gsMapData<T> & map_ori,
                         const gsMapData<T> & map_def,
                               gsMaterialData<T> & data) const
    {
        if (map_def.flags==0)
            _computeStrains_impl(map_ori,data);
        else
            _computeStrains_impl(map_ori,map_def,data);
    }

    void _computeStrains_impl(const gsMapData<T> & map_ori,
                                    gsMaterialData<T> & data) const
    {
        GISMO_ASSERT(map_ori.flags & NEED_VALUE,"The undeformed map data should contain the point values");
        GISMO_ASSERT(map_ori.flags & NEED_JACOBIAN,"The undeformed map data should contain the Jacobian matrix");
        data.dim = map_ori.points.rows();
        data.size = map_ori.points.cols();
        data.patch = (map_ori.patchId==-1) ? 0 : map_ori.patchId;
        data.m_F.resize(0,data.size);
        data.m_E.resize(0,data.size);
        data.m_parmat.resize(m_pars.size(),data.size);

        data.m_parmat.setZero();
        gsMatrix<T> pars;
        for (size_t v=0; v!=m_pars.size(); v++)
        {
            m_pars[v]->piece(data.patch).eval_into(map_ori.values[0], pars);
            data.m_parmat.row(v) = pars;
        }
    }

    void _computeStrains_impl(const gsMapData<T> & map_ori,
                              const gsMapData<T> & map_def,
                                    gsMaterialData<T> & data) const
    {
        GISMO_ASSERT(map_ori.flags & NEED_VALUE,"The undeformed map data should contain the point values");
        GISMO_ASSERT(map_ori.flags & NEED_JACOBIAN,"The undeformed map data should contain the Jacobian matrix");
        GISMO_ASSERT(map_def.flags & NEED_JACOBIAN,"The deformed map data should contain the Jacobian matrix");
        GISMO_ASSERT(map_ori.points.cols() == map_def.points.cols(),"The number of points in the undeformed and deformed maps should be the same");
        GISMO_ASSERT(map_ori.patchId == map_def.patchId,"The undeformed and deformed maps should be defined on the same patch");
        data.dim = map_ori.points.rows();
        data.size = map_ori.points.cols();
        data.patch = (map_ori.patchId==-1) ? 0 : map_ori.patchId;
        data.m_F.resize(data.dim*data.dim,data.size);
        data.m_E.resize(data.dim*data.dim,data.size);
        data.m_parmat.resize(m_pars.size(),data.size);

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
            gsAsMatrix<T, Dynamic,Dynamic> Emat = data.m_E.reshapeCol(k,data.dim,data.dim);
            Emat = 0.5 * (F.transpose() * F - I);
        }

        for (size_t v=0; v!=m_pars.size(); v++)
        {
            m_pars[v]->piece(data.patch).eval_into(map_ori.values[0], pars);
            data.m_parmat.row(v) = pars;
        }
    }

    void _computeGeometricData(const gsFunctionSet<T> * undeformed,
                               const gsFunctionSet<T> * deformed,
                               const index_t patch,
                               const gsMatrix<T>& u,
                                     gsMapData<T> & map_ori,
                                     gsMapData<T> & map_def) const
    {
        if (deformed==nullptr || deformed->nPieces()==0)
            _computeGeometricData_impl(undeformed,patch,u,map_ori,map_def);
        else
            _computeGeometricData_impl(undeformed,deformed,patch,u,map_ori,map_def);
    }

    void _computeGeometricData_impl(const gsFunctionSet<T> * undeformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                          gsMapData<T> & map_ori,
                                          gsMapData<T> & map_def) const
    {
        map_ori.clear();
        map_def.clear();
        GISMO_ASSERT(undeformed->domainDim()==u.rows(),"Geometric dimension and the point dimension are not the same!\n"<<undeformed->domainDim()<<"!="<<u.rows());
        map_ori.flags = NEED_VALUE;
        map_ori.points = u;
        map_ori.patchId = patch;
        undeformed->function(patch).computeMap(map_ori);
    }

    void _computeGeometricData_impl(const gsFunctionSet<T> * undeformed,
                                    const gsFunctionSet<T> * deformed,
                                    const index_t patch,
                                    const gsMatrix<T>& u,
                                          gsMapData<T> & map_ori,
                                          gsMapData<T> & map_def) const
    {
        map_ori.clear();
        map_def.clear();
        GISMO_ASSERT(undeformed->domainDim()==deformed->domainDim(),"Geometric dimension and the point dimension are not the same!\n"<<undeformed->domainDim()<<"!="<<deformed->domainDim());
        GISMO_ASSERT(undeformed->domainDim()==u.rows(),"Geometric dimension and the point dimension are not the same!\n"<<undeformed->domainDim()<<"!="<<u.rows());

        map_ori.flags = NEED_VALUE;
        map_ori.flags|= NEED_JACOBIAN;
        map_def.flags = NEED_JACOBIAN;
        map_def.points = map_ori.points = u;
        map_ori.patchId = map_def.patchId = patch;
        undeformed->function(patch).computeMap(map_ori);
        deformed->function(patch).computeMap(map_def);
    }

protected:

    std::vector< function_ptr > m_pars;
    function_ptr m_density;

    gsOptionList m_options;

    // Geometric data point
    // mutable util::gsThreaded< gsMaterialData<T> > m_data;

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

    typename gsMaterialBase<T>::function_ptr m_undeformed;
    typename gsMaterialBase<T>::function_ptr m_deformed;

    mutable gsMatrix<T> m_parmat;
    mutable gsMatrix<T> m_rhoMat;
    // mutable gsMatrix<T> m_jac_ori, m_jac_def;
    mutable gsMatrix<T> m_F, m_E;

    mutable short_t dim;
    mutable index_t size;
    mutable index_t patch;
};

}
