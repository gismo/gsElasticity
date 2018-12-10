/** @file gsCouplingUtilities.h

   @brief Provides functions for FSI coupling within MOTOR project,
   i.e. defining special function types for boundary conditions,
   prtojecting data between piecewise and smooth boundaries and so on.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsBoundary.h>

#include <fstream>

namespace gismo
{

template <class T>
class gsFluxFunction : public gsFunction<T>
{
public:

    gsFluxFunction(int sourcePatch,boundary::side sourceSide,gsMultiPatch<T> const & sourceGeo,
                   gsMultiPatch<T> const & sourceSolution, T alpha = 1.);

    int domainDim() const
    {
        return m_geo.domainDim();
    }

    int targetDim() const
    {
        return m_sol.targetDim();
    }

    void eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const;

protected:

    int const m_patch;
    boundary::side m_side;
    gsMultiPatch<T> const & m_geo;
    gsMultiPatch<T> const & m_sol;
    T m_alpha;

}; // class definition ends

template <class T>
class gsFluxFunctionExternal : public gsFunction<T>
{
public:

    gsFluxFunctionExternal(gsMatrix<T> const & points, gsMatrix<T> const & normGrads, T conductivity = 1.);

    int domainDim() const
    {
        return m_points.rows();
    }

    int targetDim() const
    {
        return 1;
    }

    void eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const;

protected:

    gsMatrix<T> const & m_points;
    gsMatrix<T> const & m_grads;
    T m_cond;

}; // class definition ends

template <class T>
class gsGradFunction : public gsFunction<T>
{
public:

    gsGradFunction(int sourcePatch,gsMultiPatch<T> const & sourceSolution);

    int domainDim() const
    {
        return m_sol.domainDim();
    }

    int targetDim() const
    {
        return m_sol.targetDim()*m_sol.domainDim();
    }

    void eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const;

protected:

    int const m_patch;
    gsMultiPatch<T> const & m_sol;

}; // class definition ends

template <class T>
class gsPhysGradFunction : public gsFunction<T>
{
public:

    gsPhysGradFunction(int sourcePatch,gsMultiPatch<T> const & sourceGeo,
                       gsMultiPatch<T> const & sourceSolution);

    int domainDim() const
    {
        return m_geo.domainDim();
    }

    int targetDim() const
    {
        return m_sol.targetDim()*m_geo.domainDim();
    }

    void eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const;

protected:

    int const m_patch;
    gsMultiPatch<T> const & m_geo;
    gsMultiPatch<T> const & m_sol;

}; // class definition ends

template <class T>
class gsDetFunction : public gsFunction<T>
{
public:

    gsDetFunction(int sourcePatch,gsMultiPatch<T> const & sourceGeo);

    int domainDim() const
    {
        return m_geo.domainDim();
    }

    int targetDim() const
    {
        return 1;
    }

    void eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const;

protected:

    int const m_patch;
    gsMultiPatch<T> const & m_geo;

}; // class definition ends

template <class T>
void gsExtractNormalGradients(gsMultiPatch<T> const & sourceGeo, gsMultiPatch<T> const & sourceSol,
                              std::vector<int> patches, std::vector<boundary::side> sides,
                              gsMatrix<T> & points, gsMatrix<T> & normGrads);
template <class T>
void gsProjectPiecewiseLin(gsMatrix<T> const & point, gsMatrix<T> const & set, T & ratio, int & index);

template <class T>
void gsDistPointSegment(gsMatrix<T> const & A, gsMatrix<T> const & B, gsMatrix<T> const & C,
                        T & dist, T & ratio);

} // namespace gismo
