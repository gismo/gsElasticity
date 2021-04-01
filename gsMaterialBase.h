/** @file gsMaterialBase.h

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

template <class T>
class gsMaterialBase
{

public:
    virtual ~gsMaterialBase() {};

    /**
     * @brief      { function_description }
     *
     * @param[in]  u       { parameter_description }
     * @param      result  The result
     */
    virtual void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result,
                                  index_t k = 0) const = 0;

    /**
     * @brief      { function_description }
     *
     * @param[in]  u       { parameter_description }
     * @param      result  The result
     */
    virtual void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result,
                                  index_t k = 0) const = 0;

};

}
