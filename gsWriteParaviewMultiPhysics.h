/** @file gsWriteParaviewMultiPhysics.h

    @brief Allows to write several fields defined on the same geometry
    in one file, making it easier to operate with them inside Paraview.
    Ideally should be a part of gismoIO module,

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy (TU Kaiserslautern)
    Inspired by gsWriteParaview.h by A. Mantzaflaris
*/
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsParaviewCollection.h>

#define NS 1000

#Hello world
#Bye

namespace gismo
{
/// \brief Write a file containing several fields defined on the same geometry to ONE paraview file
///
/// \param fields a map of field pointers
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
template<class T>
void gsWriteParaviewMultiPhysics(std::map<std::string, const gsField<T> *> fields, std::string const & fn,
                     unsigned npts=NS, bool mesh = false, bool ctrlNet = false);


/// \brief Write a file containing several fields defined on the same geometry to ONE paraview file
/// and adds it as a timestep to a Paraview collection
/// \param fields a map of field pointers
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
template<class T>
void gsWriteParaviewMultiPhysicsTimeStep(std::map<std::string, const gsField<T> *> fields, std::string const & fn,
                                         gsParaviewCollection & collection, int time, unsigned npts=NS);


/// \brief Extract and evaluate geometry and the fields for a single patch
///
/// \param fields a map of field pointers
/// \param patchNum a number of patch
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
template<class T>
void gsWriteParaviewMultiPhysicsSinglePatch(std::map<std::string, const gsField<T> *> fields,
                                const unsigned patchNum,
                                std::string const & fn,
                                unsigned npts);


/// \brief Utility function to actually write prepaired matrices with data into Paraview file
///
/// \param points a matrix with space-points to plot
/// \param data a map of matrices with field evaluations to plotfilename where paraview file is written
/// \param np a vector containg the data range info
/// \param fn filename where paraview file is written
template<class T>
void gsWriteParaviewMultiTPgrid(gsMatrix<T> const& points,
                                std::map<std::string, gsMatrix<T> >& data,
                                const gsVector<index_t> & np,
                                std::string const & fn);

}

#undef NS

