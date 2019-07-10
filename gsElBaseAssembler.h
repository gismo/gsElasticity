/** @file gsElBaseAssembler.h

    @brief Base class for assemblers of gsElasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

namespace gismo
{

/** @brief Extends the gsAssembler class by adding functionality necessary for a general nonlinear solver.
 * Potentially, can be merged back into gsAssembler.
 */
template <class T>
class gsElBaseAssembler : public gsAssembler<T>
{
public:
    virtual void assemble(const gsMatrix<T> & solutionVector) = 0;

    static gsOptionList defaultOptions()
    {
        gsOptionList opt = gsAssembler<T>::defaultOptions();
        opt.addReal("ForceScaling","Force scaling parameter",1.);
        opt.addReal("DirichletScaling","Dirichlet BC scaling parameter",1.);
        return opt;
    }

};

} // namespace ends
