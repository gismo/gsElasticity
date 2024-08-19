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

template <class T, enum gsMaterialOutput out> 
class gsMaterialEvalSingle;

template <class T, enum gsMaterialOutput out>
class gsMaterialEval : public gsFunctionSet<T>
{

public:

    /// Constructor
    gsMaterialEval( const gsMaterialContainer<T> & materialMatrices,
                    const gsFunctionSet<T>       * deformed)
    :
    m_materials(materialMatrices),
    m_deformed(deformed)
    {
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            GISMO_ASSERT(materialMatrices.piece(p)->initialized(),"Material matrix "<<p<<" is incomplete!");

        this->_makePieces();
    }

    /// Constructor
    gsMaterialEval(       gsMaterialBase<T> * materialMatrix,
                    const gsFunctionSet<T>    * deformed)
    :
    m_materials(deformed->nPieces()),
    m_deformed(deformed)
    {
        // GISMO_ASSERT(materialMatrix->initialized(),"Material matrix is incomplete!");
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            m_materials.set(p,materialMatrix);
        this->_makePieces();
    }

    /// Constructor
    gsMaterialEval( const gsMaterialContainer<T> & materialMatrices,
                    const gsFunctionSet<T>       * undeformed,
                    const gsFunctionSet<T>       * deformed)
    :
    m_materials(materialMatrices),
    m_deformed(deformed)
    {
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            GISMO_ASSERT(materialMatrices.piece(p)->initialized(),"Material matrix "<<p<<" is incomplete!");

        this->_makePieces(undeformed);
    }

    /// Constructor
    gsMaterialEval(       gsMaterialBase<T> * materialMatrix,
                    const gsFunctionSet<T>    * undeformed,
                    const gsFunctionSet<T>    * deformed)
    :
    m_materials(deformed->nPieces()),
    m_deformed(deformed)
    {
        // GISMO_ASSERT(materialMatrix->initialized(),"Material matrix is incomplete!");
        for (index_t p = 0; p!=deformed->nPieces(); ++p)
            m_materials.set(p,materialMatrix);
        this->_makePieces(undeformed);
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

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }

protected:
    void _makePieces()
    {
        m_pieces.resize(m_deformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialEvalSingle<T,out>(p,m_materials.piece(p),m_deformed);
    }

    void _makePieces(const gsFunctionSet<T> * undeformed)
    {
        m_pieces.resize(m_deformed->nPieces());
        for (size_t p = 0; p!=m_pieces.size(); ++p)
            m_pieces.at(p) = new gsMaterialEvalSingle<T,out>(p,m_materials.piece(p),undeformed,m_deformed);
    }

protected:
    gsMaterialContainer<T> m_materials;
    const gsFunctionSet<T> * m_deformed;
    gsMatrix<T> m_z;
    mutable std::vector<gsMaterialEvalSingle<T,out> *> m_pieces;
}; //gsMaterialEval

/**
 * @brief      This class serves as the evaluator of materials, based on \ref gsMaterialBase
 *
 * @tparam     T     Real type
 * @tparam     out   Output type (see \ref MaterialOutput)
 *
 * @ingroup    Elasticity
 */
template <class T, enum gsMaterialOutput out>
class gsMaterialEvalSingle : public gsFunction<T>
{
public:

    /// Constructor
    gsMaterialEvalSingle(   index_t patch,
                            gsMaterialBase<T>     * materialMatrix,
                            const gsFunctionSet<T>  * deformed)
    :
    m_pIndex(patch),
    m_material(materialMatrix)
    {
        m_material->setDeformed(deformed);
    }

    /// Constructor
    gsMaterialEvalSingle(   index_t patch,
                            gsMaterialBase<T>     * materialMatrix,
                            const gsFunctionSet<T>  * undeformed,
                            const gsFunctionSet<T>  * deformed)
    :
    m_pIndex(patch),
    m_material(materialMatrix)
    {
        m_material->setDeformed(undeformed);
        m_material->setDeformed(deformed);
    }

    short_t domainDim() const {return m_material->dim();}

    /**
     * @brief      Target dimension
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const {return targetDim_impl<out>();}


private:

    /// Implementation of \ref targetDim for density (TODO), energy
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Psi, short_t>::type targetDim_impl() const 
    { 
        return 1; 
    };

    /// Implementation of \ref targetDim for strain, stress
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::E ||
                            _out==gsMaterialOutput::S  , short_t>::type targetDim_impl() const 
    { 
        const short_t d = m_material->dim();
        return d*d; 
    };

    /// Implementation of \ref targetDim for strain, stress
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::E_voigt ||
                            _out==gsMaterialOutput::S_voigt  , short_t>::type targetDim_impl() const 
    { 
        const short_t d = m_material->dim();
        return d*(d+1)/2;
    };

    /// Implementation of \ref targetDim for matrix C (size = ((d+1)*d/2) x ((d+1)*d/2))
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::C  , short_t>::type targetDim_impl() const 
    { 
        return math::pow((m_material->dim()+1)*m_material->dim()/2,2); 
    };

protected:
    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

public:

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        this->eval_into_impl<out>(u,result);
    }

private:
    /// Specialisation of \ref eval_into for density (TODO), energy
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Psi   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_energy_into(m_pIndex,u,result);
    }

    /// Specialisation of \ref eval_into for the strain tensor
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::E     , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_strain_into(m_pIndex,u,result);
    }

    /// Specialisation of \ref eval_into for the strain tensor
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::E_voigt   , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_material->eval_strain_into(m_pIndex,u,tmp); //tmp gives strain tensor (dxd)
        const short_t d = m_material->dim();
        result.resize(d*(d+1)/2,u.cols());
        
        for (index_t k = 0; k < u.cols(); k++)
        {
            // gsAsMatrix<T,Dynamic,Dynamic> E = tmp.reshapeCol(k,d,d); // in tensor notation
            gsMatrix<T> E = tmp.reshapeCol(k,d,d); // in tensor notation
            gsVector<T> E_voigt; // voigt strain
            voigtStress(E_voigt,E);
            result.col(k) = E_voigt;
        }

    }

    /// Specialisation of \ref eval_into for the stress tensor
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::S     , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_stress_into(m_pIndex,u,result);
    }

    /// Specialisation of \ref eval_into for the material tensor
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::C     , void>::type eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_material->eval_matrix_into(m_pIndex,u,result);
    }

protected:
    index_t m_pIndex;
    gsMaterialBase<T> * m_material;
};

} // namespace
