/** @file gsBasePde.h

    @brief IMHO, a useless class, but it is necessary to use the gsAssembler class.
    Contains proper information for elasticity and Navier-Stokes solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/


#pragma once

#include <gsPde/gsPde.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{

template<class T>
class gsBasePde : public gsPde<T>
{

public:

    /// Default constructor
    gsBasePde( ) { }

    /// Constructor
    gsBasePde(const gsMultiPatch<T>         &domain,
              const gsBoundaryConditions<T> &bc,
              const gsPiecewiseFunction<T>  &rhs)
    : gsPde<T>(domain,bc), m_rhs(rhs)
    {
        m_unknownDim.setOnes(m_rhs.targetDim());
    }

    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const { return 1; }
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const { return os; }
    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }
    virtual int numUnknowns() const     {return m_rhs.targetDim();}

protected:
    using gsPde<T>::m_unknownDim;
    gsPiecewiseFunction<T> m_rhs;
};

} // namespace gismo
