/** @file gsMaterialEval.h

    @brief Visitor class for volumetric integration of the linear elasticity system.

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

namespace gismo
{

#include <gsCore/gsFunction.h>

enum class gsMaterialOutput : short_t
{
    Vector   = 0,
    Matrix   = 1,
};

template <class T, enum gsMaterialOutput out>
class gsMaterialEval : public gsFunction<T>
{

public:
    /// Constructor
    gsMaterialEval(   gsMaterialBase<T> * materialMatrix,
                            short_t dim)
    :
    m_materialMat(materialMatrix),
    m_dim(dim)
    {
        m_pIndex = 0;
    }

    /// Domain dimension
    short_t domainDim() const {return m_dim;}

    /**
     * @brief      Target dimension
     *
     * For a scalar (e.g. density) the target dimension is 1, for a vector (e.g. stress tensor in Voight notation) the target dimension is 3 and for a matrix (e.g. the material matrix) the target dimension is 9, which can be reshaped to a 3x3 matrix.
     *
     * @return     Returns the target dimension depending on the specified type (scalar, vector, matrix etc.)
     */
    short_t targetDim() const { return targetDim_impl<out>(); }

    /// Implementation of piece, see \ref gsFunction
    const gsFunction<T> & piece(const index_t p) const
    {
        m_piece = new gsMaterialEval(*this);
        m_piece->setPatch(p);
        return *m_piece;
    }

    /// Sets the patch index
    void setPatch(index_t p) {m_pIndex = p; }

    /// Destructor
    ~gsMaterialEval() { delete m_piece; }

    /// Implementation of eval_into, see \ref gsFunction
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        eval_into_impl<out>(u,result);
    }

private:
    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Vector, short_t>::type
    targetDim_impl() const
    {
        return m_dim;
    };

    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Matrix, short_t>::type
    targetDim_impl() const
    {
        return m_dim*m_dim;
    };

    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out!=gsMaterialOutput::Matrix &&
                            _out!=gsMaterialOutput::Vector, short_t>::type
    targetDim_impl() const
    {
        GISMO_NO_IMPLEMENTATION
    };

    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Vector, void>::type
    eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_materialMat->eval_vector_into(u,result,m_pIndex);
    };

    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out==gsMaterialOutput::Matrix, void>::type
    eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_materialMat->eval_matrix_into(u,result,m_pIndex);
    };

    template<enum gsMaterialOutput _out>
    typename std::enable_if<_out!=gsMaterialOutput::Matrix &&
                            _out!=gsMaterialOutput::Vector, void>::type
    eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION
    };

protected:
    gsMaterialBase<T> * m_materialMat;
    mutable gsMaterialEval<T,out> * m_piece;
    index_t m_pIndex;
    size_t m_dim;
};

}
