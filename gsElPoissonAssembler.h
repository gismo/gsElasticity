/** @file gsElPoissonAssembler.h

    @brief Provides stiffness matrix for Poisson's equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>

namespace gismo
{

template <class T>
class gsElPoissonAssembler : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;

    gsElPoissonAssembler(const gsMultiPatch<T> & patches,
                         const gsMultiBasis<T> & basis,
                         const gsBoundaryConditions<T> & bconditions,
                         const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the mass matrix
    virtual void assemble(bool saveEliminationMatrix = false);

    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDDoFs) {assemble();}

    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & displacement) const;

protected:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_options;
    using Base::m_system;
    using Base::m_ddof;
    using Base::eliminationMatrix;
};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElPoissonAssembler.hpp)
#endif
