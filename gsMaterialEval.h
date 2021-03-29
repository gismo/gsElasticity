/** @file gsMaterialEval.h

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

#include <gsCore/gsFunction.h>
#include <gsElasticity/gsMaterialBase.h>

namespace gismo
{

template <class T, bool matrix = true>
class gsMaterialEval : public gsFunction<T>
{

public:
    gsMaterialEval( gsMaterialBase<T> * materialMatrix)
    :
    m_materialMat(materialMatrix),
    m_piece(NULL)
    {
        m_pIndex = 0;
    }

    short_t domainDim() const {return m_materialMat->dim();}
    short_t targetDim() const {return targetDim_impl<matrix>();}

private:
    template <bool _matrix>
    typename std::enable_if<_matrix, short_t>::type
    targetDim_impl() const
    {
        return m_materialMat->size()*m_materialMat->size();
    }

    template <bool _matrix>
    typename std::enable_if<!_matrix, short_t>::type
    targetDim_impl() const
    {
        return m_materialMat->size();
    }

public:
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
        eval_into_impl<matrix>(u,result);
    }

private:
    template <bool _matrix>
    typename std::enable_if<_matrix, void>::type
    eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_materialMat->eval_matrix_into(u,result);
    }

    template <bool _matrix>
    typename std::enable_if<!_matrix, void>::type
    eval_into_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_materialMat->eval_vector_into(u,result);
    }

protected:
    gsMaterialBase<T> * m_materialMat;
    mutable gsMaterialEval<T> * m_piece;
    index_t m_pIndex;
};

}