/** @file gsALE.h

    @brief Implementation of mesh deformation method for partitioned FSI solver.

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
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsBaseAssembler.h>

namespace gismo
{

template <class T>
class gsBaseAssembler;
template <class T>
class gsIterative;

template <class T>
class gsALE
{
public:
    gsALE(gsMultiPatch<T> & geometry, const gsMultiPatch<T> & displacement,
          const gsBoundaryInterface & interfaceStr2Mesh, ale_method::method method);

    /// default option list. used for initialization
    static gsOptionList defaultOptions();

    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }

    /// number of degrees of freedom
    index_t numDofs() const {return assembler->numDofs();}

    /// construct ALE displacement field
    void constructSolution(gsMultiPatch<T> & solution) const;

    /// update mesh to comply with the current displacement field
    index_t updateMesh();

    /// save module state
    void saveState();

    /// recover module state from saved state
    void recoverState();

    /// get FSI interface container to access patch sides
    const gsBoundaryInterface & interfaceStr2Mesh() { return interface;}

protected:
    void initialize();

    /// update mesh using HE, LE or BHE methods
    index_t linearMethod();

    /// update mesh using IHE, ILE or IBHE methods
    index_t linearIncrementalMethod();

    /// update mesh using TINE or TINE_StVK methods
    index_t nonlinearMethod();

protected:
    /// outer displacement field that drives the mesh deformation
    const gsMultiPatch<T> & disp;
    /// mapping between patch sides of the fluid and solid
    const gsBoundaryInterface & interface;
    /// mesh deformation method
    ale_method::method methodALE;
    /// option list
    gsOptionList m_options;
    /// assembler
    typename gsBaseAssembler<T>::uPtr assembler;
    /// nonlinear solver
    typename gsIterative<T>::uPtr solverNL;
    /// current ALE displacement field
    gsMultiPatch<T> ALEdisp;
    /// initialization flag
    bool initialized;


    /// saved state
    bool hasSavedState;
    gsMultiPatch<T> ALEdispSaved;
};

} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALE.hpp)
#endif
