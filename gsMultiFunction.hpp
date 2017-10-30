#pragma once

namespace gismo
{

template <class T>
gsMultiFunction<T>::gsMultiFunction()
    : m_dim(-1),
      m_parDim(-1)
{

}

//template <class T>
//gsMultiFunction<T>::~gsMultiFunction()
//{
    //freeAll(m_functions);
//}

template <class T>
//void gsMultiFunction<T>::addFunction(gsFunction<T> * f)
void gsMultiFunction<T>::addFunction(std::unique_ptr<gsFunction<T> > f)
{
    if ( m_dim == -1 )
    {
        m_dim = f->domainDim();
        m_parDim = f->targetDim();
    } else
    {
        GISMO_ASSERT( m_dim == f->domainDim() && m_parDim == f->targetDim() ,
                      "Tried to add a function of different dimension in a multifunction." );
    }
    m_functions.push_back(std::move(f));
}

} //namespace ends
