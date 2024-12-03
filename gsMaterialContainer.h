/** @file gsMaterialContainer.h

    @brief Provides a container for material matrices for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsElasticity/gsMaterialBase.h>

namespace gismo
{


/**
 * @brief      This class serves as the evaluator of materials, based on \ref gsMaterialBase
 *
 * @tparam     T     Real type
 * @tparam     out   Output type (see \ref gsMaterialOutput)
 *
 * @ingroup    Elasticity
 */
template <class T>
class gsMaterialContainer // change name to PtrContainer
{
public:
    typedef typename std::vector<typename gsMaterialBase<T>::Ptr> Container;

    typedef typename Container::iterator       iterator;
    typedef typename Container::const_iterator const_iterator;

    /// Shared pointer for gsMaterialContainer
    typedef memory::shared_ptr< gsMaterialContainer > Ptr;

    /// Unique pointer for gsMaterialContainer
    typedef memory::unique_ptr< gsMaterialContainer > uPtr;
public:

    /// Constructor
    gsMaterialContainer( index_t size = 0 )
    {
        m_container.resize(size);
        // To do: initialize with null pointers
    }

    /// Constructor
    gsMaterialContainer( const gsMaterialBase<T> * mat, index_t size = 1 )
    {
        m_container.reserve(size);
        for (index_t i=0; i!=size; i++)
            this->add(mat);
    }

    gsMaterialContainer(const gsMaterialContainer & other)
    {
        // for (index_t k=0; k!=other.m_container.size(); k++)
        //     add(memory::make_unique(other.m_container.at(k)));
        m_container = give(other.m_container);
    }

    ~gsMaterialContainer()
    {
        // freeAll(m_container);
    }

    /// Add a material matrix by copying argument
    void add(const gsMaterialBase<T> & mat)
    {
        m_container.push_back( memory::make_shared(mat.clone().release()) );
    }

    ///\brief Add a material matrix from a gsMaterialBase<T>::uPtr
    void add(const gsMaterialBase<T> * mat)
    {
        m_container.push_back( memory::make_shared_not_owned(mat) );
    }

    /// Set a material matrix by copying argument
    void set(const index_t i, const gsMaterialBase<T> & mat)
    {
        m_container[i] = memory::make_shared(mat.clone().release());
    }

    ///\brief Set a material matrix from a gsMaterialBase<T>::uPtr
    void set(const index_t i, const gsMaterialBase<T> * mat)
    {
        m_container[i] = memory::make_shared_not_owned(mat);
    }

    ///\brief Set a material matrix from a gsMaterialBase<T>::uPtr
    void set(const index_t i, const typename gsMaterialBase<T>::Ptr mat)
    {
        m_container[i] = mat;
    }

    gsMaterialBase<T> * piece(const index_t k) const
    {
        return m_container.at(k).get();
    }

    index_t size() const {return m_container.size();}

    std::ostream &print(std::ostream &os) const
    {
        os << "Material container with "<<m_container.size() <<" pieces.\n";
        return os;
    }

    friend std::ostream & operator<<(std::ostream & os, const gsMaterialContainer & mc)
    {
        return mc.print(os);
    }

    /// Clear all function pointers
    void clear()
    {
        m_container.clear();
    }

protected:
    Container m_container;

};

// TODO: Check gsKLShell
// namespace internal
// {

// /// @brief get a MaterialContainer from XML data
// ///
// /// \ingroup KLShell
// template<class T>
// class gsXml< gsMaterialContainer<T> >
// {
// private:
//     gsXml() { }
//     typedef gsMaterialContainer<T> Object;

// public:
//     GSXML_COMMON_FUNCTIONS(gsMaterialContainer<T>);
//     static std::string tag ()  { return "MaterialContainer"; }
//     static std::string type () { return ""; }

//     GSXML_GET_POINTER(Object);

//     static void get_into(gsXmlNode * node,Object & obj)
//     {
//         const int size = atoi(node->first_attribute("size")->value());

//         // Read material inventory
//         int count = countByTag("Material", node);
//         std::vector<typename gsMaterialBase<T>::Ptr> mat(count);
//         for (gsXmlNode * child = node->first_node("Material"); child; child =
//                 child->next_sibling("Material"))
//         {
//             const int i = atoi(child->first_attribute("index")->value());
//             mat[i] = memory::make_shared(gsXml<gsMaterialBase<T>>::get(child));
//         }

//         obj = gsMaterialContainer<T>(size);
//         for (gsXmlNode * child = node->first_node("group"); child;
//                 child = child->next_sibling("group"))
//         {
//             const int mIndex = atoi(child->first_attribute("material")->value());
//             std::istringstream group_str;
//             group_str.str( child->value() );

//             for(int patch; ( gsGetInt(group_str,patch)); )
//                 obj.set(patch,mat[mIndex]);
//         }

//     }

//     static gsXmlNode * put (const Object & obj,
//                             gsXmlTree & data)
//     {
//         GISMO_ERROR("Writing gsMaterialContainer to Xml is not implemented");
//         // gsWarn<<"Writing gsMaterialContainer to Xml is not implemented\n";
//         // gsXmlNode * result;
//         // return result;
//         // return putMaterial< Object >( obj,data );
//     }
// };

// } // namespace internal

} // namespace gismo
