/** @file gsPhaseFieldAssembler.h

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

enum PForder
{
    Second = 2,
    Fourth = 4
};

enum PFmode
{
    AT1 = 0,
    AT2 = 1
};


namespace gismo
{

template <class T>
class gsPhaseFieldAssemblerBase
{
public:

    virtual ~gsPhaseFieldAssemblerBase() = default;

    virtual void initialize() = 0;

    virtual void setSpaceBasis(const gsFunctionSet<T> & spaceBasis) = 0;

    /**
     * @brief Assembles the matrix \Phi (see Greco et al 2024)
     */
    virtual void assembleMatrix() = 0;

    /**
     * @brief Assembles the residual vector
     */
    virtual void assembleResidual(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC) = 0;


    /**
     * @brief Assembles int_\Omega \Psi()*w*w.tr()*meas(G)
     */
    virtual void assemblePsiMatrix(const gsFunctionSet<T> & PSI) = 0;

    /**
     * @brief Assembles int_\Omega \Psi()*w*meas(G)
     */
    virtual void assemblePsiVector(const gsFunctionSet<T> & PSI) = 0;

    virtual const gsSparseMatrix<T> & matrix() const = 0;

    virtual void matrix_into(gsSparseMatrix<T> & out) = 0;

    virtual const gsMatrix<T>       & rhs()    const = 0;

    virtual void rhs_into(gsMatrix<T> & out) = 0;

    virtual index_t numDofs() const = 0;

    virtual gsOptionList & options() = 0;

    virtual void setOptions(gsOptionList & options) = 0;

    virtual void constructSolution(gsMatrix<T>     & Cvec,
                                   gsMultiPatch<T> & C) const = 0;

    virtual void constructSolution(gsMatrix<T>         & Cvec,
                                      gsMappedSpline<2,T> & C) const = 0;

    virtual void constructSolution(const gsMultiPatch<T> & C,
                                      gsMatrix<T>     & Cvec) const = 0;

};

template <class T, enum PForder order, enum PFmode mode>
class gsPhaseFieldAssembler : public gsPhaseFieldAssemblerBase<T>
{
public:

/**
 * @brief      Constructs a new biharmonic assembler
 *
 * @param[in]  mp     A multi-patch
 * @param[in]  mb     A multi-basis
 * @param[in]  force  The force
 * @param[in]  bcs    The bcs
 */
gsPhaseFieldAssembler(const gsMultiPatch<T> & mp,
                      const gsMultiBasis<T> & mb,
                      const gsBoundaryConditions<T> & bcs
                      );

/// Default empty constructor
gsPhaseFieldAssembler() { }

/// Copy constructor (makes deep copy)
gsPhaseFieldAssembler( const gsPhaseFieldAssembler& other )
{
    operator=(other);
}

/// Move constructor
gsPhaseFieldAssembler( gsPhaseFieldAssembler&& other )
{
    operator=(give(other));
}

/// Assignment operator
gsPhaseFieldAssembler& operator= ( const gsPhaseFieldAssembler& other );

/// Move assignment operator
gsPhaseFieldAssembler& operator= ( gsPhaseFieldAssembler&& other );

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

    // /**
    //  * @brief Assembles the mass matrix separately
    //  */
    // void assembleMassMatrix();

    /**
     * @brief Assembles the matrix \Phi (see Greco et al 2024)
     */
    void assembleMatrix() override;

    /**
     * @brief Assembles the residual
     * @param C The solution
     * @param DC The time-derivative of the solution
     */
    void assembleResidual(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC) override;

    /**
     * @brief Assembles int_\Omega \Psi()*w*w.tr()*meas(G)
     */
    void assemblePsiMatrix(const gsFunctionSet<T> & PSI) override;

    /**
     * @brief Assembles int_\Omega \Psi()*w*meas(G)
     */
    void assemblePsiVector(const gsFunctionSet<T> & PSI) override;

    // /**
    //  * @brief Assembles the Nitsche vector for boundary conditions separately
    //  * @param C The solution
    //  * @param DC The time-derivative of the solution
    //  */
    // void assembleNitscheVector(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC);

    // /**
    //  * @brief Assembles the Nitsche matrix for boundary conditions separately
    //  * @note  This term does not depend on the solution, hence it can be assembled once
    //  */
    // void assembleNitscheMatrix();

    /**
     * @brief Returns a handle to the latest assembled matrix
     * @return The matrix
     */
    const gsSparseMatrix<T> & matrix() const override { return m_assembler.matrix();  }

    /**
     * @brief Moves the latest assembled matrix from the assembler to the output
     * @param out The output matrix
     */
    void matrix_into(gsSparseMatrix<T> & out) override { m_assembler.matrix_into(out); }

    /**
     * @brief Returns a handle to the latest assembled right-hand side
     * @return The right-hand side
     */
    const gsMatrix<T>       & rhs()    const  override{ return m_assembler.rhs();     }

    /**
     * @brief Moves the latest assembled right-hand side from the assembler to the output
     * @param out The output right-hand side
     */
    void rhs_into(gsMatrix<T> & out) override { m_assembler.rhs_into(out); }

    /**
     * @brief Returns the number of degrees of freedom
     * @return The number of degrees of freedom
     */
    index_t numDofs() const override { return m_assembler.numDofs(); };

    /**
     * @brief Returns a handle to the options stored in the class
     * @return The options
     */
    gsOptionList & options() override {return m_options;}

    /**
     *  @brief Set the options from an option list. Ignores unknown options
     */
    void setOptions(gsOptionList & options) override;

    /**
     * @brief Constructs a multi-patch solution from a solution vector
     * @param Cvec The solution vector
     * @param C The solution
     */
    void constructSolution(gsMatrix<T>     & Cvec,
                           gsMultiPatch<T> & C) const override;

    /**
     * @brief Constructs a spline solution from a solution vector
     * @param Cvec The solution vector
     * @param C The solution
     */
    void constructSolution(gsMatrix<T>         & Cvec,
                           gsMappedSpline<2,T> & C) const override;

    /**
     * @brief Constructs a solution vector from a multi-patch solution
     * @param C The solution
     * @param Cvec The solution vector
     */
    void constructSolution(const gsMultiPatch<T> & C,
                                 gsMatrix<T>     & Cvec) const override;

private:

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT1), void>::type
    _assembleMatrix_impl();

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT1), void>::type
    _assembleMatrix_impl();

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT2), void>::type
    _assembleMatrix_impl();

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT2), void>::type
    _assembleMatrix_impl();



    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT1), void>::type
    assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC);

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT1), void>::type
    assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC);

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT2), void>::type
    assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC);

    template <enum PForder _order, enum PFmode _mode>
    typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT2), void>::type
    assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC);

protected:

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;
    typedef typename gsExprAssembler<T>::element     element;

    mutable index_t m_continuity;
    mutable T m_l0, m_Gc;
    // mutable T m_penalty;

    gsExprAssembler<T> m_assembler;
    gsExprEvaluator<T> m_evaluator;

    gsMultiPatch<T>           m_patches;
    mutable gsMultiBasis<T>   m_basis;
    const gsFunctionSet<T> *  m_spaceBasis;
    gsBoundaryConditions<T>   m_bcs;
    bool m_initialized;

    mutable gsOptionList m_options;

  }; // class gsPhaseFieldAssembler

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsPhaseFieldAssembler
   */
  void pybind11_init_gsPhaseFieldAssembler(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPhaseFieldAssembler.hpp)
#endif
