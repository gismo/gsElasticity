/** @file gsMaterialEval.h

    @brief Evaluator for `gsMaterialBase` objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

#include <gsElasticity/gsMaterialUtils.h>
#include <gsElasticity/gsMaterialContainer.h>
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsCore/gsFunction.h>

namespace gismo
{

template <class T, enum gsMaterialOutput out, bool voigt>
class gsMaterialEvalSingle;

template <class T, enum gsMaterialOutput out, bool voigt = true>
class gsMaterialEval : public gsFunctionSet<T>
{
    typedef typename gsMaterialBase<T>::function_ptr function_ptr;

public:

    /// Constructor
    gsMaterialEval( const gsMaterialContainer<T> & materialMatrices,
                    const gsFunctionSet<T>       & undeformed)
    :
    m_materials(materialMatrices),
    m_undeformed(undeformed),
    m_deformed(nullptr)
    {
        this->_makePieces(m_undeformed);
    }

    /// Constructor
    gsMaterialEval(       gsMaterialBase<T> * materialMatrix,
                    const gsFunctionSet<T>    & undeformed)
    :
    m_materials(undeformed.nPieces()),
    m_undeformed(undeformed),
    m_deformed(nullptr)
    {
        for (index_t p = 0; p!=undeformed.nPieces(); ++p)
            m_materials.set(p,materialMatrix);
        this->_makePieces(m_undeformed);
    }

    /// Constructor
    gsMaterialEval( const gsMaterialContainer<T> & materialMatrices,
                    const gsFunctionSet<T>       & undeformed,
                    const gsFunctionSet<T>       & deformed)
    :
    m_materials(materialMatrices),
    m_undeformed(undeformed),
    m_deformed(deformed)
    {
        this->_makePieces(m_undeformed,m_deformed);
    }

    /// Constructor
    gsMaterialEval(       gsMaterialBase<T> * materialMatrix,
                    const gsFunctionSet<T>    & undeformed,
                    const gsFunctionSet<T>    & deformed)
    :
    m_materials(deformed.nPieces()),
    m_undeformed(memory::make_shared_not_owned(undeformed.clone().release())),
    m_deformed(memory::make_shared_not_owned(deformed.clone().release()))
    {
        for (index_t p = 0; p!=deformed.nPieces(); ++p)
            m_materials.set(p,materialMatrix);
        this->_makePieces(m_undeformed,m_deformed);
    }

    /// Destructor
    ~gsMaterialEval()
    {
        freeAll(m_pieces);
    }

    /// Domain dimension
    short_t domainDim() const {return this->piece(0).domainDim();}

    /**
     * @brief      Target dimension
     *
     * For a scalar (e.g. density) the target dimension is 1, for a vector (e.g. stress tensor in Voight notation) the target dimension is 3 and for a matrix (e.g. the material matrix) the target dimension is 9, which can be reshaped to a 3x3 matrix.
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const { return this->piece(0).targetDim(); }

    /// Implementation of piece, see \ref gsFunction
    const gsFunction<T> & piece(const index_t p) const
    {
        return *m_pieces[p];
    }

    /// Implementation of nPieces(), see \ref gsFunctionSet
    index_t nPieces() const override{return m_pieces.size();}

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

protected:
    void _makePieces(function_ptr undeformed)
    {
        m_pieces.resize(undeformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialEvalSingle<T,out,voigt>(p,m_materials.piece(p),undeformed,nullptr);
    }

    void _makePieces(function_ptr undeformed, function_ptr deformed)
    {
        m_pieces.resize(deformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialEvalSingle<T,out,voigt>(p,m_materials.piece(p),undeformed,deformed);
    }

protected:
    gsMaterialContainer<T> m_materials;
    function_ptr m_undeformed;
    function_ptr m_deformed;
    mutable std::vector<gsMaterialEvalSingle<T,out,voigt> *> m_pieces;
}; //gsMaterialEval

/**
 * @brief      This class serves as the evaluator of materials, based on \ref gsMaterialBase
 *
 * @tparam     T     Real type
 * @tparam     out   Output type (see \ref MaterialOutput)
 * @tparam     voigt   Voigt notation (true: Voigt, false: tensor)
 *
 * @ingroup    Elasticity
 */
template <class T, enum gsMaterialOutput out, bool voigt = true>
class gsMaterialEvalSingle : public gsFunction<T>
{
    typedef typename gsMaterialBase<T>::function_ptr function_ptr;

public:

    /// Constructor
    gsMaterialEvalSingle(   index_t patch,
                            gsMaterialBase<T>   * materialMatrix,
                            function_ptr undeformed,
                            function_ptr deformed)
    :
    m_pIndex(patch),
    m_material(materialMatrix),
    m_undeformed(undeformed),
    m_deformed(deformed),
    m_dim(m_undeformed->domainDim())
    {
    }

    short_t domainDim() const {return m_dim;}

    /**
     * @brief      Target dimension
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const {return targetDim_impl<out,voigt>();}


private:

    /// Implementation of \ref targetDim for density (TODO), energy
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::Psi, short_t>::type targetDim_impl() const
    {
        return 1;
    };

    /// Implementation of \ref targetDim for strain, stress
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<(_out==gsMaterialOutput::E ||
                            _out==gsMaterialOutput::S) && !_voigt , short_t>::type targetDim_impl() const
    {
        return m_dim*m_dim;
    };

    /// Implementation of \ref targetDim for strain, stress
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<(_out==gsMaterialOutput::E ||
                            _out==gsMaterialOutput::S) && _voigt , short_t>::type targetDim_impl() const
    {
        return m_dim*(m_dim+1)/2;
    };

    /// Implementation of \ref targetDim for matrix C (size = ((d+1)*d/2) x ((d+1)*d/2)) -- in Voigt
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C  && _voigt, short_t>::type targetDim_impl() const
    {
        return math::pow((m_dim+1)*m_dim/2,2);
    };

    /// Implementation of \ref targetDim for matrix C (size = ((d+1)*d/2) x ((d+1)*d/2))
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C  && !_voigt, short_t>::type targetDim_impl() const
    {
        GISMO_NO_IMPLEMENTATION;
    };

        /// Implementation of \ref targetDim for matrix C (size = ((d+1)*d/2) x ((d+1)*d/2)) -- in Voigt
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C_pos  && _voigt, short_t>::type targetDim_impl() const
    {
        return math::pow((m_dim+1)*m_dim/2,2);
    };

    /// Implementation of \ref targetDim for matrix C (size = ((d+1)*d/2) x ((d+1)*d/2))
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C_pos  && !_voigt, short_t>::type targetDim_impl() const
    {
        GISMO_NO_IMPLEMENTATION;
    };

protected:
    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

public:

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        this->eval_into_impl<out, voigt>(u,result);
    }

private:
    /// Specialisation of \ref eval_into for density (TODO), energy
    template<enum gsMaterialOutput _out, bool>
    typename std::enable_if<_out==gsMaterialOutput::Psi, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_energy_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
    }

    /// Specialisation of \ref eval_into for the strain tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::E  && !_voigt, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_strain_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
    }

    /// Specialisation of \ref eval_into for the strain tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::E && _voigt   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_material->eval_strain_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
        result.resize(m_dim*(m_dim+1)/2,u.cols());
        calculate_voigt_strain(tmp, m_dim, result);
        // for (index_t k = 0; k < u.cols(); k++)
        // {
        //     // gsAsMatrix<T,Dynamic,Dynamic> E = tmp.reshapeCol(k,d,d); // in tensor notation
        //     gsMatrix<T> E = tmp.reshapeCol(k,d,d); // in tensor notation
        //     gsVector<T> E_voigt; // voigt strain
        //     voigtStress(E_voigt,E);
        //     result.col(k) = E_voigt;
        // }

    }

    /// Specialisation of \ref eval_into for the stress tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::S && !_voigt, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_stress_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
    }


    /// Specialisation of \ref eval_into for the strain tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::S && _voigt, void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_material->eval_stress_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
        result.resize(m_dim*(m_dim+1)/2,u.cols());
        calculate_voigt_stress(tmp, m_dim, result);
    }

    /// Specialisation of \ref eval_into for the material tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C && _voigt , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_matrix_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
    }

    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C && !_voigt , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// Specialisation of \ref eval_into for the material tensor
    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C_pos && _voigt , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_matrix_pos_into(*m_undeformed,*m_deformed,m_pIndex,u,result);
    }

    template<enum gsMaterialOutput _out, bool _voigt>
    typename std::enable_if<_out==gsMaterialOutput::C_pos && !_voigt , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

public:


public:

    /// Precompute the geometric data
    void precompute(const gsMatrix<T>& u)
    {
        m_material->precomputeData(m_pIndex,u);
    }

    /// Get the geometric data

protected:
    index_t m_pIndex;
    gsMaterialBase<T> * m_material;
    function_ptr m_undeformed;
    function_ptr m_deformed;
    const short_t m_dim;
};

} // namespace
