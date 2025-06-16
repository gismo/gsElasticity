/** @file gsSolidAssembler.h

    @brief Provides assembler for a (planar) Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------

    Author(s): H.M. Verhelst, L. Venta Viñuela
*/

#pragma once

#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{

template <class T>
class gsSolidAssemblerBase
{
public:

    virtual ~gsSolidAssemblerBase() = default;

    virtual void initialize() = 0;

    virtual void setSpaceBasis(const gsFunctionSet<T> & spaceBasis) = 0;

    virtual void assemble() = 0;

    virtual void assemble(const gsMatrix<T> & uvec) = 0;

    virtual void assembleMass() = 0;

    virtual const gsSparseMatrix<T> matrix() const = 0;

    virtual void matrix_into(gsSparseMatrix<T> & out) = 0;

    virtual const gsMatrix<T> rhs() const = 0;

    virtual void rhs_into(gsMatrix<T> & out) = 0;

    virtual index_t numDofs() const = 0;

    virtual gsOptionList & options() = 0;

    virtual void setOptions(gsOptionList & options) = 0;

    virtual void constructSolution(gsMatrix<T> & uvec,
                                   gsMultiPatch<T> & displacement) const = 0;

    void constructDisplacement(gsMatrix<T> & uvec,
                               gsMultiPatch<T> & displacement) const
    {
        this->constructSolution(uvec, displacement);
    }

    virtual void constructDeformed(gsMatrix<T> & uvec,
                                   gsMultiPatch<T> & deformed) const = 0;

    virtual void constructSolution(const gsMultiPatch<T> & displacement, gsMatrix<T> & uvec) const = 0;


};

template <short_t DIM, class T, class Material>
class gsSolidAssembler : public gsSolidAssemblerBase<T>
{
public:

gsSolidAssembler(const gsMultiPatch<T> & mp,
                 const gsMultiBasis<T> & mb,
                 const gsBoundaryConditions<T> & bcs,
                 const Material *         material,
                 const gsFunctionSet<T> * g=nullptr);

/// Default empty constructor
gsSolidAssembler() { }

/// Copy constructor (makes deep copy)
gsSolidAssembler( const gsSolidAssembler& other )
{
    operator=(other);
}

/// Move constructor
gsSolidAssembler( gsSolidAssembler&& other )
{
    operator=(give(other));
}

/// Assignment operator
gsSolidAssembler& operator= ( const gsSolidAssembler& other );

/// Move assignment operator
gsSolidAssembler& operator= ( gsSolidAssembler&& other );

protected:

    void _defaultOptions();

    void _getOptions();

public:

    void initialize() override;

    /**
     * @brief Overwrites the basis to be used as a space,
     *        and keeps the basis defined in the constructor
     *        for integration
     * @note  This option is usually called for assembly on \ref gsMappedBasis
     * @param spaceBasis The basis to be used for the space
     */
    void setSpaceBasis(const gsFunctionSet<T> & spaceBasis) override;

    void assemble() override;

    void assemble(const gsMatrix<T> & uvec) override;

    void assembleMass() override;

    const gsSparseMatrix<T> matrix() const override { return m_assembler.matrix(); }

    void matrix_into(gsSparseMatrix<T> & out) override { m_assembler.matrix_into(out); }

    const gsMatrix<T> rhs() const override { return m_assembler.rhs(); }

    void rhs_into(gsMatrix<T> & out) override { m_assembler.rhs_into(out); }

    index_t numDofs() const override { return m_assembler.numDofs(); };

    gsOptionList & options() override {return m_options;}

    void setOptions(gsOptionList & options) override;

    virtual void constructSolution(gsMatrix<T> & uvec, gsMultiPatch<T> & displacement) const;

    virtual void constructDeformed(gsMatrix<T> & uvec, gsMultiPatch<T> & deformed) const;

    virtual void constructSolution(const gsMultiPatch<T> & displacement, gsMatrix<T>     & uvec) const;

protected:

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;
    typedef typename gsExprAssembler<T>::variable    variable;
    typedef typename gsExprAssembler<T>::element     element;

    mutable index_t m_continuity;
    mutable bool m_smallStrain;

    gsExprAssembler<T> m_assembler;
    gsExprEvaluator<T> m_evaluator;

    gsMultiPatch<T>                 m_patches;
    mutable gsMultiBasis<T>         m_basis;
    const gsFunctionSet<T> *        m_spaceBasis;
          gsBoundaryConditions<T>   m_bcs;
    const gsFunctionSet<T> *        m_forcing;
    const Material *                m_material;
    bool  m_initialized;

    mutable gsOptionList m_options;

  }; // class gsSolidAssembler

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsSolidAssembler
   */
  void pybind11_init_gsSolidAssembler(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSolidAssembler.hpp)
#endif
