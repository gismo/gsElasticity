/** @file gsElPoissonPde.h

   @brief Allows to introduce a scaling constant in front of the stiffness matrix.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsPde/gsPoissonPde.h>

namespace gismo
{

template <class T>
class gsElPoissonPde : public gsPoissonPde<T>
{
public:

    gsElPoissonPde(const gsMultiPatch<T>         &domain,
                   const gsBoundaryConditions<T> &bc,
                   const gsPiecewiseFunction<T>  &rhs,
                         T                       k)
        : gsPoissonPde<T>(domain,bc,rhs,NULL), m_k(k)
    {

    }

    virtual T k() const { return m_k; }
protected:
    T m_k;

}; // class definition ends

} // namespace ends

