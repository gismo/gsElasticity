/** @file gsWriteParaviewMultiPhysics.cpp

    @brief Provides implementation for gsWriteParaviewMultiPhysics.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy (TU Kaiserslautern)
*/

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsIO/gsWriteParaview.hpp>

#define PLOT_PRECISION 5

namespace gismo
{


template<class T>
void gsWriteParaviewMultiPhysics(std::map<std::string, const gsField<T>*> fields,
                                 std::string const & fn,
                                 unsigned npts, bool mesh)
{
    const unsigned n = fields.begin()->second->nPatches();
    gsParaviewCollection collection(fn);
    std::string fileName = fn.substr(fn.find_last_of("/\\")+1); // file name without a path

    for ( unsigned i=0; i < n; ++i )
    {
        const gsBasis<> & dom = fields.begin()->second->isParametrized() ?
            fields.begin()->second->igaFunction(i).basis() : fields.begin()->second->patch(i).basis();

        writeSinglePatchMultiField( fields, i, fn + util::to_string(i), npts);
        collection.addPart(fileName + util::to_string(i), ".vts");

        if ( mesh )
        {
            writeSingleCompMesh(dom, fields.begin()->second->patch(i), fn + util::to_string(i) + "_mesh");
            collection.addPart(fileName + util::to_string(i) + "_mesh", ".vtp");
        }
    }
    collection.save();
}

template<class T>
void writeSinglePatchMultiField(std::map<std::string,const gsField<T> *> fields,
                                const unsigned patchNum,
                                std::string const & fn,
                                unsigned npts)
{
    const gsFunction<> & geometry = fields.begin()->second->patch(patchNum);
    const unsigned n = geometry.targetDim();
    const unsigned d = geometry.domainDim();

    gsMatrix<> ab = geometry.support();
    gsVector<> a = ab.col(0);
    gsVector<> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a,b,npts);
    gsMatrix<> pts = gsPointGrid(a,b,np);

    gsMatrix<> eval_geo = geometry.eval(pts);
    std::map<std::string, gsMatrix<> > data;
    for (typename std::map<std::string,const gsField<T> *>::iterator it = fields.begin(); it != fields.end(); it++)
    {
        data[it->first] = it->second->isParametric() ?
                    it->second->function(patchNum).eval(pts) : it->second->function(patchNum).eval(eval_geo);
        if ( data[it->first].rows() > 1 )
        {
            data[it->first].conservativeResize(3,eval_geo.cols() );
            data[it->first].bottomRows( 3-d ).setZero();
        }
    }

    if (3 -d > 0)
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }
    else if (d > 3)
    {
        gsWarn<< "Cannot plot 4D data.\n";
        return;
    }

    if ( 3 - n > 0 )
    {
        eval_geo.conservativeResize(3,eval_geo.cols() );
        eval_geo.bottomRows(3-n).setZero();
    }
    else if (n > 3)
    {
        gsWarn<< "Data is more than 3 dimensions.\n";
    }

    for (typename std::map<std::string, gsMatrix<> >::iterator it = data.begin(); it != data.end(); it++)
    {
        if ( it->second.rows() > 1 )
        {
            it->second.conservativeResize(3,eval_geo.cols() );
            it->second.bottomRows( 3-d ).setZero();
        }
    }

    gsWriteParaviewMultiTPgrid(eval_geo, data, np.template cast<index_t>(), fn);
}

template<class T>
void gsWriteParaviewMultiTPgrid(gsMatrix<T> const& points,
                                std::map<std::string, gsMatrix<T> >& data,
                                const gsVector<index_t> & np,
                                std::string const & fn)
{
    const int n = points.rows();

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "
         << (np.size()>2 ? np(2)-1 : 0) <<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "
         << (np.size()>2 ? np(2)-1 : 0) <<"\">\n";

    file <<"<PointData>\n";
    for (typename std::map<std::string, gsMatrix<T> >::iterator it = data.begin(); it != data.end(); it++)
    {
        file <<"<DataArray type=\"Float32\" Name=\""<< it->first <<"\" format=\"ascii\" NumberOfComponents=\""<< ( it->second.rows()==1 ? 1 : 3) <<"\">\n";
        if ( it->second.rows()==1 )
            for ( index_t j=0; j<it->second.cols(); ++j)
                file<< it->second.at(j) <<" ";
        else
        {
            for ( index_t j=0; j<it->second.cols(); ++j)
            {
                for ( index_t i=0; i!=it->second.rows(); ++i)
                    file<< it->second(i,j) <<" ";
                for ( index_t i=it->second.rows(); i<3; ++i)
                    file<<"0 ";
            }
        }
        file <<"</DataArray>\n";
    }
    file <<"</PointData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    for ( index_t j=0; j<points.cols(); ++j)
    {
        for ( index_t i=0; i!=n; ++i)
            file<< points(i,j) <<" ";
        for ( index_t i=n; i<3; ++i)
            file<<"0 ";
    }
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";

    file.close();
}

TEMPLATE_INST
void gsWriteParaviewMultiPhysics(std::map<std::string, const gsField<real_t>* > fields, std::string const & fn,
                     unsigned npts, bool mesh);

TEMPLATE_INST
void writeSinglePatchMultiField(std::map<std::string,const gsField<real_t>* > fields,
                                const unsigned patchNum,
                                std::string const & fn,
                                unsigned npts);

TEMPLATE_INST
void gsWriteParaviewMultiTPgrid(gsMatrix<real_t> const& points,
                                std::map<std::string, gsMatrix<real_t> >& data,
                                const gsVector<index_t> & np,
                                std::string const & fn);


}

#undef PLOT_PRECISION
