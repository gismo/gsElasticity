/** @file gsBaseAssembler.h

    @brief Base class for assemblers of gsElasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
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
class gsBaseAssembler : public gsAssembler<T>
{
public:
    typedef memory::shared_ptr<gsBaseAssembler> Ptr;
    typedef memory::unique_ptr<gsBaseAssembler> uPtr;

    /// Assembles the tangential linear system for Newton's method given the current solution
    /// in the form of free and fixed/Dirichelt degrees of freedom.
    /// Checks if the current solution is valid (Newton's solver can exit safely if invalid).
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDDoFs) = 0;

    /// assembly procedure for linear problems
    virtual void assemble(bool saveEliminationMatrix = false) {};

    /// Returns number of free degrees of freedom
    virtual int numDofs() const { return gsAssembler<T>::numDofs(); }

    /// Constructs solution as a gsMultiPatch object from the solution vector and fixed DoFs
    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDDofs,
                                   gsMultiPatch<T> & result,
                                   const gsVector<index_t> & unknowns) const;

    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDDofs,
                                   gsMultiPatch<T> & result) const {};

    //--------------------- DIRICHLET BC SHENANIGANS ----------------------------------//

    /** @brief Set Dirichet degrees of freedom on a given side of a given patch from a given matrix.
     *
     * A usual source of degrees of freedom is another geometry where it is known that the corresponding
     * bases match. The function is not safe in that it assumes a correct numbering of DoFs in a given
     * matrix. To allocate space for these DoFs in the assembler, add an empty/zero Dirichlet boundary condition
     * to gsBoundaryCondtions container that is passed to the assembler constructor.
     */
    virtual void setFixedDofs(index_t patch, boxSide side, const gsMatrix<T> & ddofs, bool oneUnk = false);

    /// set all fixed degrees of freedom
    virtual void setFixedDofs(const std::vector<gsMatrix<T> > & ddofs);

    /// get fixed degrees of freedom corresponding to a given part of the bdry.
    /// each column of the resulting matrix correspond to one variable/component of the vector-valued vairable
    virtual void getFixedDofs(index_t patch, boxSide side, gsMatrix<T> & ddofs) const;

    /// get the size of the Dirichlet vector for elimination
    virtual index_t numFixedDofs() const;

    /// @brief Eliminates new Dirichelt degrees of fredom
    virtual void eliminateFixedDofs();

    //virtual void modifyDirichletDofs(size_t patch, boxSide side, const gsMatrix<T> & ddofs);

    //--------------------- OTHER ----------------------------------//

    virtual void setRHS(const gsMatrix<T> & rhs) {m_system.rhs() = rhs;}

    virtual void setMatrix(const gsSparseMatrix<T> & mat) {m_system.matrix() = mat;}

protected:
    using gsAssembler<T>::m_pde_ptr;
    using gsAssembler<T>::m_bases;
    using gsAssembler<T>::m_system;
    using gsAssembler<T>::m_ddof;

    gsSparseMatrix<T> eliminationMatrix;
    gsMatrix<T> rhsWithZeroDDofs;
};

} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBaseAssembler.hpp)
#endif
