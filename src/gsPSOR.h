/** @file gsPSOR.h

    @brief Solves a Symmetric Linear Complementary Problem using the
           Projected Successive Over-Relaxation (PSOR) method.

           Solves
              \f[
                \begin{aligned}
                    x^T (M x + q) &= 0,\\
                    M x + q &>= 0, \\
                    x &>= 0.
                \end{aligned}
                \f]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. M. Verhelst
*/
#pragma once

namespace gismo
{

template<class T = real_t>
class gsPSOR
{
public:
    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    explicit gsPSOR( const gsSparseMatrix<T> & mat)
    : m_mat(mat),
      m_max_iters(1000),
      m_tolU(1e-10),
      m_tolNeg(1e-10),
      m_tolPos(1e-10),
      m_num_iter(-1),
      m_rhs_norm(-1),
      m_error{{-1,-1,-1}},
      m_verbose(false)
    {
        GISMO_ASSERT(m_mat.rows()     == m_mat.cols(),     "The matrix is not square."                     );
        this->defaultOptions();
    }

    /// @brief Returns a list of default options
    void defaultOptions()
    {
        m_options.addInt   ("MaxIterations"    , "Maximum number of iterations", 1000       );
        m_options.addReal  ("tolU"  , "Tolerance for the error criteria on the relative residual error",      1e-10      );
        m_options.addReal  ("tolNeg", "Tolerance for the error criteria on the relative residual error",      1e-10      );
        m_options.addReal  ("tolPos", "Tolerance for the error criteria on the relative residual error",      1e-10      );
        m_options.addSwitch("Verbose", "Print information about the iteration", false);
    }

    const gsOptionList & options() const { return m_options; } ///< Access the options
          gsOptionList & options()       { return m_options; } ///< Access the options

    /// @brief Set the options based on a gsOptionList
    void setOptions(const gsOptionList & options)
    {
        m_options.update(options,gsOptionList::ignoreIfUnknown);
    }

    // Register options
    void getOptions()
    {
        m_max_iters        = m_options.askInt   ("MaxIterations"  , m_max_iters        );
        m_tolU             = m_options.askReal  ("tolU"           , m_tolU             );
        m_tolNeg           = m_options.askReal  ("tolNeg"         , m_tolNeg           );
        m_tolPos           = m_options.askReal  ("tolPos"         , m_tolPos           );
        m_verbose          = m_options.askSwitch("Verbose"        , m_verbose          );
    }

    void solve(const gsMatrix<T> & rhs, gsMatrix<T> & x)
    {
        this->getOptions();
        if (initIteration(rhs, x)) return;

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;
            if (step(rhs,x)) break;
        }
        GISMO_ASSERT(m_num_iter < m_max_iters, "The PSOR solver did not converge. It = " << m_num_iter << "/" << m_max_iters << ", errors = " << m_error[0] << ", " << m_error[1] << ", " << m_error[2]);
        finalizeIteration(x);
    }

    bool initIteration(const gsMatrix<T> & rhs, gsMatrix<T> & x)
    {
        GISMO_ASSERT( rhs.cols() == 1,
                      "Iterative solvers only work for single column right-hand side." );
        GISMO_ASSERT( rhs.rows() == m_mat.rows(),
                      "The right-hand side does not match the matrix: "
                      << rhs.rows() <<"!="<< m_mat.rows() );

        m_num_iter = 0;

        m_rhs_norm = rhs.norm();

        if (0 == m_rhs_norm) // special case of zero rhs
        {
            x.setZero(rhs.rows(),rhs.cols()); // for sure zero is a solution
            m_error = {0.,0.,0.};
            return true; // iteration is finished
        }

        if ( 0 == x.size() ) // if no initial solution, start with zeros
            x.setZero(rhs.rows(), rhs.cols());
        else
        {
            GISMO_ASSERT( x.cols() == 1,
                      "Iterative solvers only work for single right-hand side and solution." );
            GISMO_ASSERT( x.rows() == m_mat.cols(),
                      "The initial guess does not match the matrix: "
                      << x.rows() <<"!="<< m_mat.cols() );
        }
        return false; // iteration is not finished
    }


    bool step( const gsMatrix<T> & rhs, gsMatrix<T> & x )
    {
        // update solution @ previous iteration x_im1
        gsMatrix<T> x_im1 = x;

        std::vector<T> err_d_vec (m_mat.cols(), 0);
        std::vector<T> err_a0_vec(m_mat.cols(), 0);
        std::vector<T> err_ap_vec(m_mat.cols(), 0);
        T err_d = 0;
        T err_a0 = 0;
        T err_ap = 0;
        T act;

        T Qx = 0;   // scalar product Q*x_im1
        T Ldx = 0;  // scalar product L*(x_i-x_im1)
        T D_m1 = 0; // scalar         Q_{jj}^-1

        index_t Irow, Jcol;

/*
        NOTE ON PARALLELIZATION:

        This algorithm cannot be parallelized since the strict upper-triangular terms
        depend on x(Irow), which might be updated by another thread.
*/
        // loop over solution components (i.e. Q matrix columns)
        for (index_t j = 0; j < m_mat.outerSize(); ++j)
        {
            Qx = Ldx = D_m1 = 0; // reset
            // loop over matrix rows
            typename gsSparseMatrix<T>::iterator mIt = m_mat.begin(j);
            Jcol = mIt.col();
            for ( ; mIt; ++mIt)
            {
                Irow = mIt.index();
                Qx += mIt.value() * x_im1(Irow);
                // Diagonal term
                if (Irow == mIt.col())
                    D_m1 = math::pow(mIt.value(),-1.0);
                // strict upper triangular term
                if (Irow < mIt.col())
                    Ldx += mIt.value() * (x(Irow) - x_im1(Irow));
            }

            x(Jcol) = math::max(0.0, x_im1(Jcol) - D_m1 * (Qx + rhs(Jcol) + Ldx));
        }

        // compute error
        for (index_t j = 0; j < m_mat.outerSize(); ++j)
        {
            Qx = 0;
            typename gsSparseMatrix<T>::iterator mIt = m_mat.begin(j);
            Jcol = mIt.col();
            for ( ; mIt; ++mIt)
            {
                Irow = mIt.index();
                Qx += mIt.value() * x(Irow);
            }
            act = Qx + rhs(Jcol);

            err_d = math::max(err_d, math::abs(x(Jcol) - x_im1(Jcol))); // non-parallel
            if ( x(Jcol) > 0)
                err_ap = math::max(err_ap, math::abs(act) ); // non-parallel
            else
                err_a0 = math::min(err_a0, act); // non-parallel
        }

        if (m_verbose) gsInfo<<"PSOR It. "<<m_num_iter<<": ||x|| = "<<x.norm()<<", err_d = "<<err_d<<", err_ap = "<<err_ap<<", err_a0 = "<<err_a0<<"\n";

        // return (err_d < m_tol);
        m_error = {err_d, err_ap, err_a0};
        return (err_d < m_tolU && err_ap < m_tolPos && err_a0 < m_tolNeg);
    }

    virtual void finalizeIteration( const gsMatrix<T> & ) {}          ///< Some post-processing might be required

    index_t numIter() const { return m_num_iter; } ///< Returns the number of iterations

    std::array<T,3> errors() const { return m_error; }

private:
    const gsSparseMatrix<T> & m_mat;
    index_t m_max_iters;
    T m_tolU, m_tolNeg, m_tolPos;
    index_t m_num_iter;
    T m_rhs_norm;
    std::array<T,3> m_error;
    bool m_verbose;
    gsOptionList m_options;
};

} // namespace gismo
