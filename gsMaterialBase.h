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
// #include <gsUtils/gsThreaded.h>

namespace gismo
{

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
    m_undeformed(nullptr),
    m_deformed(nullptr),
    m_density(nullptr)
    {
    }

    gsMaterialBase( const gsFunctionSet<T> & mp,
                    const gsFunctionSet<T> & mp_def,
                    const gsFunctionSet<T> & Density)
:
    m_undeformed(memory::make_shared(mp.clone().release())),
    m_deformed(memory::make_shared(mp_def.clone().release())),
    m_density(memory::make_shared(Density.clone().release()))
    {
    }

    gsMaterialBase( const gsFunctionSet<T> * mp,
                    const gsFunctionSet<T> * mp_def,
                    const gsFunctionSet<T> * Density)
    :
    m_undeformed(memory::make_shared_not_owned(mp)),
    m_deformed(memory::make_shared_not_owned(mp_def)),
    m_density(memory::make_shared_not_owned(Density))
    {
    }

    gsMaterialBase( const function_ptr & mp,
                    const function_ptr & mp_def,
                    const function_ptr & Density)
    :
    m_undeformed(mp),
    m_deformed(mp_def),
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
            m_undeformed = other.m_undeformed;
            m_deformed = other.m_deformed;
            m_options = other.m_options;

            m_pars = other.m_pars;
            m_density = other.m_density;
        }
        return *this;
    }

    gsMaterialBase<T> & operator= ( gsMaterialBase<T> && other )
    {
        m_undeformed = give(other.m_undeformed);
        m_deformed = give(other.m_deformed);
        m_options = give(other.m_options);


        m_pars = give(other.m_pars);
        m_density = give(other.m_density);
        return *this;
    }

    /**
     * @brief      Sets the default options
     */
    virtual inline void defaultOptions()
    { GISMO_NO_IMPLEMENTATION; }

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
    { GISMO_NO_IMPLEMENTATION; }

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
    virtual inline void eval_deformation_gradient_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }


    /**
     * @brief      Evaluates the strain tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return     
     */
    virtual inline void eval_strain_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }


    /**
     * @brief      Evaluates the stress tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return     
     */
    virtual inline void eval_stress_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the material tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return     
     */
    virtual inline void eval_matrix_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the undamaged positive material tensor at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return     
     */
    virtual inline void eval_matrix_pos_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }

    /**
     * @brief      Evaluates the energy at a set of given points
     *
     * @param[in]  patch  The patch
     * @param[in]  u      The in-plane shell coordinates to be eveluated on
     *
     * @return     
     */
    virtual inline void eval_energy_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T> & result) const
    { GISMO_NO_IMPLEMENTATION; }

    // //////// CONTINUE HERE
    // void gradient_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Fres, index_t k = 0) const
    // {
    //     gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
    //     gsMapData<T> map_def = _getMapData(deformed,u,k);

    //     gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
    //     gsMatrix<T> physDispJac;
    //     Fres.resize(u.rows()*u.rows(),u.cols());
    //     for (index_t i=0; i!=u.cols(); i++)
    //     {
    //         physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());
    //         // deformation gradient F = I + du/dx
    //         gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,u.rows(),u.rows());
    //         F = I + physDispJac;
    //     }
    // }

    // void gradient_into(const gsMatrix<T>& u, gsMatrix<T>& Fres, index_t k = 0) const
    // {
    //     this->gradient_into(m_deformed,u,Fres,k);
    // }

    // void strain_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Eres, index_t k = 0) const
    // {
    //     gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
    //     gsMapData<T> map_def = _getMapData(deformed,u,k);

    //     gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
    //     gsMatrix<T> physDispJac, F, RCG;
    //     Eres.resize(u.rows()*u.rows(),u.cols());
    //     for (index_t i=0; i!=u.cols(); i++)
    //     {
    //         physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());

    //         // deformation gradient F = I + du/dx = dx/dX
    //         F = I + physDispJac;
    //         // Right Cauchy Green strain, C = F'*F
    //         RCG = F.transpose() * F;
    //         // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
    //         gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
    //         E = 0.5 * (RCG - I);
    //     }
    // }

    // void strain_into(const gsMatrix<T>& u, gsMatrix<T>& Eres, index_t k = 0) const
    // {
    //     this->strain_into(m_deformed,u,Eres,k);
    // }

    // void allStrains_into(const gsFunctionSet<T> * deformed, const gsMatrix<T>& u, gsMatrix<T>& Fres, gsMatrix<T>& Eres, index_t k = 0) const
    // {
    //     gsMapData<T> map_ori = _getMapData(m_undeformed,u,k);
    //     gsMapData<T> map_def = _getMapData(deformed,u,k);

    //     gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
    //     gsMatrix<T> physDispJac, F, RCG, E;
    //     Fres.resize(u.rows()*u.rows(),u.cols());
    //     Eres.resize(u.rows()*u.rows(),u.cols());
    //     for (index_t i=0; i!=u.cols(); i++)
    //     {
    //         physDispJac = map_def.jacobian(i)*(map_ori.jacobian(i).cramerInverse());
    //         // deformation gradient F = I + du/dx = dx/dX
    //         gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,u.rows(),u.rows());
    //         F = I + physDispJac;
    //         // Right Cauchy Green strain, C = F'*F
    //         RCG = F.transpose() * F;
    //         // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
    //         gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
    //         E = 0.5 * (RCG - I);
    //     }
    // }

    // void allStrains_into(const gsMatrix<T>& u, gsMatrix<T>& Fres, gsMatrix<T>& Eres, index_t k = 0) const
    // {
    //     this->allStrains_into(m_deformed,u,Fres,Eres,k);
    // }

    /// Sets the density
    virtual inline void setDensity(function_ptr Density) { m_density = Density; }
    /// Sets the density
    virtual void setDensity(const gsFunctionSet<T> & Density)
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

    virtual void setUndeformed(const gsFunctionSet<T> * undeformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(undeformed);
        m_undeformed = f_ptr;
    }
    virtual void setDeformed(const gsFunctionSet<T> * deformed)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(deformed);
        m_deformed = f_ptr;
    }

    virtual void setUndeformed(const function_ptr undeformed) {m_undeformed = undeformed; }
    virtual void setDeformed(const function_ptr deformed) {m_deformed = deformed; }

    const function_ptr getUndeformed() const { return m_undeformed; }
    const function_ptr getDeformed()  const { return m_deformed; }

    virtual bool hasUndeformed() const
    { return m_deformed!=nullptr; }
    virtual bool hasDeformed() const { return m_deformed!=nullptr; }

    virtual bool initialized() const
    { GISMO_NO_IMPLEMENTATION; }

protected:

    function_ptr m_undeformed;
    function_ptr m_deformed;

    std::vector< function_ptr > m_pars;
    function_ptr m_density;

    gsOptionList m_options;

public:

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //must be present whenever the class contains fixed size matrices
#   undef Eigen

};
}
