/** @file gsPartitionedFSI.h

    @brief Partitioned FSI solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once
#include <gsIO/gsOptionList.h>
#include <gsElasticity/gsBaseUtils.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

template <class T>
class gsPartitionedFSI
{
public:

};

} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPartitionedFSI.hpp)
#endif
