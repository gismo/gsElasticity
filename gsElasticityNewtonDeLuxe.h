/** @file gsElasticityNewtonDeLuxe.h

    @brief A class providing Newton's method for nonlinear elasticity in mixed formulation.
    Supports incremental loading (ILS = incremental loading step).
    Can save every intermediate displacement field for further visualization.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template <class T>
class gsElasticityAssembler;

template <class T>
class gsElasticityNewtonDeLuxe
{
public:
    /// constructor without an initial guess
    gsElasticityNewtonDeLuxe(gsElasticityAssembler<T> & elasticityAssembler);

    static gsOptionList defaultOptions();
    gsOptionList & options() { return m_options; }
    void solve();

    /// return the last computed displacement field
    const gsMultiPatch<T> & displacement() const
    {
        GISMO_ENSURE(!displacements.empty(),"No solution computed!");
        return displacements.back();
    }
    /// return the last computed pressure field
    const gsMultiPatch<T> & pressure() const
    {
        GISMO_ENSURE(!pressures.empty(),"No solution computed!");
        return pressures.back();
    }
    /// return all saved displacement fields
    const std::vector<gsMultiPatch<T> > & allDisplacements() & { return displacements; }
    /// return all saved displacement fields
    const std::vector<gsMultiPatch<T> > & allPressures() & { return pressures; }

protected:

    void computeUpdate(bool initUpdate);

    void printStatus();

    void saveSolution(const gsMultiPatch<T> & incDisplacement,
                      const gsMultiPatch<T> & incPressure);

    void bijectivityCheck(const gsMultiPatch<T> & incDisplacement);

protected:
    /// assembler object that generates the linear system
    gsElasticityAssembler<T> & assembler;
    /// container with displacement fields; can store either only the final solution or all intermediate states as well
    std::vector<gsMultiPatch<T> > displacements;
    /// container with pressure fields; can store either only the final solution or all intermediate states as well
    std::vector<gsMultiPatch<T> > pressures;

    /// status variables
    index_t numIterations; /// number of Newton's iterations performed at the current ILS
    index_t incStep; /// current incremental loading step
    bool converged;  /// convergence status at the current ILS
    T residualNorm;
    T initResidualNorm; /// residual norm at the beginning of the current ILS
    T updateNorm;
    T initUpdateNorm; /// update vector norm at the beginning of the current ILS
    bool bijective;


    gsOptionList m_options;
};

struct newtonVerbosity2
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2  /// output everything
    };
};

struct newtonSave2
{
    enum save
    {
        onlyFinal = 0,  /// save only the final solution
        firstAndLastPerIncStep = 1,  /// save only the first and the last displacement fields at each incremental loading step
        all = 2  /// save every intermediate displacement field
    };
};

} // namespace gismo
