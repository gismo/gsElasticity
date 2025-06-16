/** @file material_expr.h

    @brief Provides an expressions for material evaluation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M.Verhelst
*/

#pragma once

#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsMaterialUtils.h>
#include <gsAssembler/gsExpressions.h>

namespace gismo
{
    namespace expr
    {
        template<class Material, bool smallStrains, enum gsMaterialOutput out, typename E>
        class material_expr  : public _expr<material_expr<Material, smallStrains, out, E> >
        {
            typedef typename Material::Scalar T;
            typename gsGeometryMap<T>::Nested_t _G;
            typename E::Nested_t _u;
            mutable gsMaterialData<T> _Mdata;
            mutable gsMatrix<T> _res;
            const short_t _d;
            const size_t _sz;
            gsMatrix<T> _I;

            std::vector<gsFeVariable<T>> _parameters;

        public:
            typedef T Scalar;
            enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

            material_expr(const gsGeometryMap<T> & G, _expr<E> const & u, std::vector<gsFeVariable<T>> & parameters)
            :
            _G(G),
            _u(u),
            _d(_G.source().domainDim()),
            _sz(_d*(_d+1)/2),
            _parameters(parameters)
            {
                _Mdata.dim = _d;
                _Mdata.size= 1;
                _Mdata.patch = 0; // not used
                _Mdata.deformationGradient.resize(_d*_d,1); // deformation gradient
                _Mdata.strain.resize(_d*_d,1); // strain tensor
                _Mdata.resizeParameters(_parameters.size()); // parameters
                _I = gsMatrix<T>::Identity(_d,_d); // identity matrix
                _Mdata.flags = NEED_MATERIAL_F | NEED_MATERIAL_E; // request deformation gradient and strain tensor
            }

            // This expression inserts the jacobian of the deformed and undeformed geometry maps into the
            // material evaluation function, computes the deformation gradient and strain tensor, and evaluates
            // the material matrix. The result is a matrix expression that can be used in the assembly process.
            const gsMatrix<T> & eval(const index_t k) const
            {
                this->_fillData(k, _Mdata); // fill the material data for the given geometry map and index k
                this-> eval_impl<out>(_Mdata,_res); // evaluate the material matrix
                return _res;
            }

            // Per basis function
            index_t rows() const { return rows_impl<out>(); }
            index_t cols() const { return cols_impl<out>(); } // dim*(dim+1)/2

            void parse(gsExprHelper<T> & evList) const
            {
                // Geometry maps
                evList.add(_G);
                _G.data().flags |= NEED_VALUE | NEED_DERIV;

                // evList.add(_u.space()); // Add the displacement expression
                // grad_expr<E>(_u).parse(evList);
                parse_impl<E>(_u,evList); // Parse _u depending on its type

                for (index_t p=0; p!=_parameters.size(); ++p)
                {
                    evList.add(_parameters[p]);
                    _parameters[p].data().flags |= NEED_VALUE;
                }
            }

            void print(std::ostream &os) const { print_impl<out>(os); }

        private:
            void _fillData(const index_t k, gsMaterialData<T> & data) const
            {
                // Compute the material data for the given geometry map and index k
                const gsMatrix<T> & jac_ori = _G.data().jacobian(k); // original jacobian
                // const gsMatrix<T> & jac_u   = _u.data().jacobian(k); // deformed jacobian

                gsMatrix<T> jac_u = gradU_impl<E>(_u, k); // evaluate deformed jacobian

                gsAsMatrix<T,Dynamic,Dynamic> deformationGradient = _Mdata.deformationGradient.reshapeCol(0, _d, _d);
                gsAsMatrix<T,Dynamic,Dynamic> strain = _Mdata.strain.reshapeCol(0, _d, _d);
                Material::template compute_strainDataFromDisplacements<smallStrains>(jac_ori, jac_u, _I, deformationGradient, strain); // compute deformation gradient and strain tensor
                for (index_t p=0; p!=_parameters.size(); ++p)
                    data.parameters[p] = _parameters[p].data().values[0].col(k); // assign parameter value to the material data
            }

            template<typename _E>
            typename std::enable_if< std::is_same<_E,gsFeSolution<Scalar>>::value , void>::type parse_impl(const _E & u, gsExprHelper<T> & evList) const
            {
                evList.add(u.space()); // Add the displacement expression
                grad_expr<E>(u).parse(evList);
            }

            template<typename _E>
            typename std::enable_if< std::is_same<_E,gsNullExpr<Scalar>>::value , void>::type parse_impl(const _E & u, gsExprHelper<T> & evList) const
            {
                //
            }

            template<typename _E>
            typename std::enable_if<!(std::is_same<_E,gsFeSolution<Scalar>>::value) &&
                                    !(std::is_same<_E,gsNullExpr<Scalar>>::value)      , void>::type parse_impl(const _E & u, gsExprHelper<T> & evList) const
            {
                GISMO_NO_IMPLEMENTATION;
            }

            template<typename _E>
            typename std::enable_if< std::is_same<_E,gsFeSolution<Scalar>>::value , gsMatrix<Scalar>>::type gradU_impl(const _E & u, const index_t k) const
            {
                grad_expr<E> jacExprU(_u); // jacobian expression for displacements
                return jacExprU.eval(k); // evaluate deformed jacobian
            }

            template<typename _E>
            typename std::enable_if< std::is_same<_E,gsNullExpr<Scalar>>::value , gsMatrix<Scalar>>::type gradU_impl(const _E & u, const index_t k) const
            {
                return gsMatrix<Scalar>::Zero(_d, _d); // return zero matrix for null expression
            }

            template<typename _E>
            typename std::enable_if<!(std::is_same<_E,gsFeSolution<Scalar>>::value) &&
                                    !(std::is_same<_E,gsNullExpr<Scalar>>::value)      , gsMatrix<Scalar>>::type gradU_impl(const _E & u, const index_t k) const
            {
                GISMO_NO_IMPLEMENTATION;
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::C , void>::type eval_impl(const gsMaterialData<T> & data, gsMatrix<T> & result) const
            {
                gsMatrix<T> tmp;
                Material::eval_matrix_into(data, tmp);
                result.resize(_sz, data.size*_sz);
                for (index_t k = 0; k < data.size; ++k)
                    result.block(0, k*_sz, _sz, _sz) = tmp.reshapeCol(k, _sz, _sz); // reshape to tensor notation
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::S , void>::type eval_impl(const gsMaterialData<T> & data, gsMatrix<T> & result) const
            {
                gsMatrix<T> tmp;
                Material::eval_stress_into(data, tmp);
                result.resize(_d, data.size*_d);
                for (index_t k = 0; k < data.size; ++k)
                    result.block(0, k*_d, _d, _d) = tmp.reshapeCol(k, _d, _d); // reshape to tensor notation
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::C , size_t>::type rows_impl() const
            {
                return _sz;
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::S , size_t>::type rows_impl() const
            {
                return _d;
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::C , size_t>::type cols_impl() const
            {
                return _sz;
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::S , size_t>::type cols_impl() const
            {
                return _d;
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::C , void>::type print_impl(std::ostream &os) const
            {
                os<< "C";
            }

            template<enum gsMaterialOutput _out>
            typename std::enable_if<_out==gsMaterialOutput::S , void>::type print_impl(std::ostream &os) const
            {
                os<< "S";
            }
        };

        template <class Material, bool smallStrains, enum gsMaterialOutput out, typename E>
        EIGEN_STRONG_INLINE
        material_expr<Material, smallStrains, out, E> const
        material(const gsGeometryMap<typename Material::Scalar> & G,
                _expr<E> const &                                  u,
                 std::vector<gsFeVariable<typename Material::Scalar>> parameters)
        {
            return material_expr<Material, smallStrains, out, E>(G, u, parameters);
        }

        template <class Material, bool smallStrains, enum gsMaterialOutput out>
        EIGEN_STRONG_INLINE
        material_expr<Material, smallStrains, out, gsNullExpr<typename Material::Scalar>> const
        material(const gsGeometryMap<typename Material::Scalar> & G,
                 std::vector<gsFeVariable<typename Material::Scalar>> parameters)
        {
            return material_expr<Material, smallStrains, out, gsNullExpr<typename Material::Scalar>>(G, gsNullExpr<typename Material::Scalar>(), parameters);
        }

    }
}
