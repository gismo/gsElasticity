/** @file gsMultiFunction.h

    @brief Provides a multi-patch container for patchwise funtions. Can be used
    for complicated functions like stresses, vorticity, etc. Allows to plot using Paraview

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/

#pragma once

#include <gsCore/gsFunctionSet.h>

namespace gismo
{

/** @brief Contains patchwise function objects derived from gsFunction.

    \tparam T coefficient type

    \ingroup Elasticity
*/
template <class T>
class gsMultiFunction : public gsFunctionSet<T>
{
public:

    gsMultiFunction();

    ~gsMultiFunction();

    int domainDim() const
    {
        return m_dim;
    }

    int targetDim() const
    {
        return m_parDim;
    }

    const gsFunctionSet<T> & piece(const index_t k) const
    {
        GISMO_ASSERT((unsigned) k < m_functions.size(), "Invalid function index requested from gsMultiFunction" );
        return *m_functions[k];
    }

    void addFunction(gsFunction<T>* f);

    void clear()
    {
        m_dim = -1;
        m_parDim = -1;
        freeAll(m_functions);
        m_functions.clear();
    }

protected:

    std::vector<gsFunction<T> *> m_functions;
    int m_dim;
    int m_parDim;

}; // gsMultiFunction ends
} // namespace ends
