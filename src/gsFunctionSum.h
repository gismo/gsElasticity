/** @file gsFunctionSum.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

#pragma once

#include<gsCore/gsFunctionSet.h>

namespace gismo
{
template<class T> class gsFunctionPieceSum;

template<class T>
class gsFunctionSum : public gsFunctionSet<T>
{
    /// Shared pointer for gsFunctionSum
    typedef memory::shared_ptr< gsFunctionSum<T> > Ptr;

    /// Unique pointer for gsFunctionExpr
    typedef memory::unique_ptr< gsFunctionSum<T> > uPtr;

public:
    gsFunctionSum()
    {}

    gsFunctionSum(const std::vector<gsFunctionSet<T> *> & functions)
    {
        for (index_t k=0; k!=m_functions.size(); k++)
            m_functions.push_back(functions[k]->clone().release());
        _init();
    }

    gsFunctionSum(const gsFunctionSet<T> * undeformed, const gsFunctionSet<T> * deformation)
    {
        m_functions.push_back(undeformed->clone().release());
        m_functions.push_back(deformation->clone().release());
        _init();
    }

    gsFunctionSum(const gsFunctionSet<T> & undeformed, const gsFunctionSet<T> & deformation)
    {
        m_functions.push_back(undeformed.clone().release());
        m_functions.push_back(deformation.clone().release());
        _init();
    }

    ~gsFunctionSum()
    {
        freeAll(m_functions);
    }

    void addFunction(typename gsFunctionSet<T>::uPtr function)
    {
        m_functions.push_back(function.release());
        _init();
    }

    void addFunction(const gsFunctionSet<T> & function)
    {
        addFunction(function.clone());
    }

    index_t nPieces() const {return m_functions[0]->nPieces();}

    short_t domainDim() const
    { return m_functions[0]->domainDim(); }

    short_t targetDim() const
    { return m_functions[0]->targetDim(); }

    const gsFunction<T> & piece(const index_t k) const override { return m_pieces[k]; }

    index_t size() const {return nPieces();}

    std::ostream &print(std::ostream &os) const
    {
        for (size_t p = 1; p!=m_pieces.size(); p++)
            gsInfo<<"Piece "<<p<<":\n"<<m_pieces[p]<<"\n";
        return os;
    }


    index_t funcSize() const {return m_functions.size();}

    typename gsFunctionSet<T>::uPtr func(const index_t i) const { return m_functions[i]->clone(); }

private:
    void _init()
    {
        GISMO_ASSERT(m_functions.size()>0,"Must give one or more function sets.");
        for (size_t p = 1; p!=m_functions.size(); p++)
        {
            GISMO_ASSERT(m_functions[p-1]->nPieces()==m_functions[p]->nPieces(),"Number of pieces does not match for function "<<p-1<<" and "<<p<<"!");
            GISMO_ASSERT(m_functions[p-1]->domainDim()==m_functions[p]->domainDim(),"Domain dimension does not match for function "<<p-1<<" and "<<p<<"!");
            GISMO_ASSERT(m_functions[p-1]->targetDim()==m_functions[p]->targetDim(),"Target dimension does not match for function "<<p-1<<" and "<<p<<"!");
        }

        for (index_t p = 0; p!=m_functions[0]->nPieces(); p++)
            m_pieces.push_back(gsFunctionPieceSum<T>(m_functions,p));
    }

protected:
    std::vector<gsFunctionSet<T> *> m_functions;
    std::vector<gsFunctionPieceSum<T> > m_pieces;

};

template<class T>
class gsFunctionPieceSum :  public gsFunction<T>
{
    /// Shared pointer for gsFunctionPieceSum
    typedef memory::shared_ptr< gsFunctionPieceSum<T> > Ptr;

    /// Auto pointer for gsFunctionExpr
    typedef memory::unique_ptr< gsFunctionPieceSum<T> > uPtr;

    using Base = gsFunction<T>;

public:
    gsFunctionPieceSum(const std::vector<gsFunctionSet<T> *> & functions,
                       const index_t index)
    :
    m_index(index),
    m_functions(functions)
    {
    }

    gsMatrix<T> support() const
    {
        gsMatrix<T> supp = m_functions[0]->piece(m_index).support();
        for (size_t p = 1; p!=m_functions.size(); p++)
            for (index_t d=0; d!=supp.rows(); d++)
            {
                GISMO_ASSERT(m_functions[p]->piece(m_index).support().rows()!=0 && m_functions[p]->piece(m_index).support().cols()!=0,"Support is empty");

                if (supp(d,0) > m_functions[p]->piece(m_index).support()(d,0))
                    supp(d,0) = m_functions[p]->piece(m_index).support()(d,0);
                if (supp(d,1) < m_functions[p]->piece(m_index).support()(d,1))
                    supp(d,1) = m_functions[p]->piece(m_index).support()(d,1);
            }
        return supp;
    }

    short_t domainDim() const
    { return m_functions[0]->domainDim(); }

    short_t targetDim() const
    {
        return m_functions[0]->targetDim();
    }



    /// Evaluates the non-zero spline functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_functions[0]->piece(m_index).eval_into(u,result);
        for (size_t p = 1; p!=m_functions.size(); p++)
        {
            m_functions[p]->piece(m_index).eval_into(u,tmp);
            result += tmp;
        }
    }

    /// Evaluates the (partial) derivatives of non-zero spline functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        gsMatrix<T> tmp;
        m_functions[0]->piece(m_index).deriv_into(u,result);
        for (size_t p = 1; p!=m_functions.size(); p++)
        {
            m_functions[p]->piece(m_index).deriv_into(u,tmp);
            result += tmp;
        }

    }

    /// Evaluates the (partial) derivatives of the nonzero spline functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        gsMatrix<T> tmp;
        m_functions[0]->piece(m_index).deriv2_into(u,result);
        for (size_t p = 1; p!=m_functions.size(); p++)
        {
            m_functions[p]->piece(m_index).deriv2_into(u,tmp);
            result += tmp;
        }

    }

    /// @brief Evaluate the nonzero spline functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, std::vector<gsMatrix<T> >& result) const
    {
        std::vector<gsMatrix<T> > tmp;
        m_functions[0]->piece(m_index).evalAllDers_into(u,n,result);
        for (size_t p = 1; p!=m_functions.size(); p++)
        {
            m_functions[p]->piece(m_index).evalAllDers_into(u,n,tmp);
            for (size_t k = 0; k!=result.size(); k++)
            {
                result[k] += tmp[k];
            }
        }
    }
    using Base::evalAllDers_into;

    GISMO_CLONE_FUNCTION(gsFunctionPieceSum)

    std::ostream &print(std::ostream &os) const
    {
        for (size_t p = 0; p!= m_functions.size(); p++)
            gsInfo<<" function "<<p<<" ("<<m_functions[p]<<"):\n"<<m_functions[p]->piece(m_index)<<"\n";
        return os;
    }

protected:
    const std::vector<gsFunctionSet<T> *> & m_functions;
    const index_t m_index;

};

} // namespace gismo
