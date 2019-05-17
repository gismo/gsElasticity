/** @file gsElasticityNewton.h

    @brief A class providing Newton's method for nonlinear elasticity.
    Supports incremental loading (ILS = incremental loading step),
    adaptive step size for bijectivity control and step size damping.
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

enum class elasticity_formulation { displacement, mixed_pressure };

// use plain enums for compatability with options
struct newton_verbosity
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2  /// output everything
    };
};

struct newton_save
{
    enum save
    {
        onlyFinal = 0,  /// save only the final solution
        firstAndLastPerIncStep = 1,  /// save only the first and the last displacement fields at each incremental loading step
        all = 2  /// save every intermediate displacement field
    };
};

template <class T>
class gsElasticityNewton
{
public:
    /// constructor without an initial guess
    gsElasticityNewton(gsElasticityAssembler<T> & elasticityAssembler,
                       elasticity_formulation form = elasticity_formulation::displacement);

    static gsOptionList defaultOptions();
    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }
    /// standard solution procedure; supports incremental loading
    void solve();
    /// adaptive solution procedure; more complex; currently only for the pure dispalcement formulation
    void solveAdaptive();

    /// return the last computed displacement field
    const gsMultiPatch<T> & displacement() const
    {
        GISMO_ENSURE(!displacements.empty(),"No solution computed!");
        return displacements.back();
    }
    /// return all saved displacement fields
    const std::vector<gsMultiPatch<T> > & allDisplacements() & { return displacements; }
    /// return the last computed pressure field
    const gsMultiPatch<T> & pressure() const
    {
        GISMO_ENSURE(formulation == elasticity_formulation::mixed_pressure, "There is no pressure in this formulation!");
        GISMO_ENSURE(!pressures.empty(),"No solution computed!");
        return pressures.back();
    }
    /// return all saved pressure fields
    const std::vector<gsMultiPatch<T> > & allPressures() & { return pressures; }
    /// use all saved displacement fields to plot the all intermediate deformed configurations of the computational domain;
    /// always plots the deformed isoparametric mesh; plots the Jacobian determinant of the deformed configuration if *numSamplingPoints* > 0
    void plotDeformation(const gsMultiPatch<T> & initDomain, std::string fileName, index_t numSamplingPoints = 10000);

protected:

    void computeUpdate(bool initUpdate);

    void printStatus();

    void saveSolution(const gsMultiPatch<T> & incDisplacement);
    void saveSolution(const gsMultiPatch<T> & incDisplacement,
                      const gsMultiPatch<T> & incPressure);

    void bijectivityCheck(const gsMultiPatch<T> & incDisplacement);

    void adaptiveHalving(gsMultiPatch<T> & incDisplacement);

    void dampedNewton(gsMultiPatch<T> & incDisplacement);

protected:
    /// assembler object that generates the linear system
    gsElasticityAssembler<T> & assembler;
    /// formulation flag for internal logic
    elasticity_formulation formulation;
    /// container with displacement fields; can store either only the final solution or all intermediate states as well
    std::vector<gsMultiPatch<T> > displacements;
    /// container with pressure fields; can store either only the final solution or all intermediate states as well
    /// only used with the mixed displacement-pressure formulation
    std::vector<gsMultiPatch<T> > pressures;

    /// status variables
    index_t numIterations; /// number of Newton's iterations performed at the current ILS
    index_t incStep; /// current incremental loading step
    bool converged;  /// convergence status at the current ILS
    T residualNorm;
    T initResidualNorm; /// residual norm at the beginning of the current ILS
    T updateNorm;
    T initUpdateNorm; /// update vector norm at the beginning of the current ILS
    index_t numAdaptHalving; /// number of adaptive stepsize halvings
    bool bijective; /// bijectivity flag for emergency stop
    T alpha; /// damping parameter

    gsOptionList m_options;
};

} // namespace gismo
