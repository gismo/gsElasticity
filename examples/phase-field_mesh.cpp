/** @file fracture_elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

    To run the script in 2D:
    ./bin/fracture_elasticity_example -f linear_elasticity_example_singlepatch_2d.xml -r 8 --plot
    ./bin/fracture_elasticity_example -f linear_elasticity_example_singlepatch_2d.xml -r 7 --plot
*/

//! [Include namespace]
#include <gismo.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsMaterialEval.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPSOR.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numHRef = 0;
    index_t numUHRef = 0;
    index_t numElev = 0;
    std::string output;
    std::string input;
    bool THB = false;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
    cmd.addInt("r", "numHRef","Number of elements in the crack size", numHRef);
    cmd.addInt("R", "numUHRef","Number of pre-refinements", numUHRef);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "input", "Input XML file", input);
    cmd.addSwitch("THB", "Use THB-splines", THB);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    GISMO_ASSERT(!input.empty(),"Input file not provided");

    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);


    gsFileData<> fd(input);
    gsInfo << "Input file "<< fd.lastPath() <<"\n";

    GISMO_ASSERT(fd.hasLabel("initial"), "Initial field not found in the input file.");
    GISMO_ASSERT(fd.hasLabel("geometry"), "Geometry not found in the input file.");
    GISMO_ASSERT(fd.hasLabel("material"), "Material parameters not found in the input file.");

    gsFunctionExpr<> initial;
    fd.getLabel("initial", initial);

    gsMultiPatch<> mp;
    fd.getLabel("geometry", mp);

    gsOptionList materialParameters;
    fd.getLabel("material", materialParameters);

    // Internal length
    real_t l0 = materialParameters.getReal("l0");
    gsMatrix<> bbox;
    mp.boundingBox(bbox);

    gsWriteParaview(mp,initial,outputdir+"initial",10000);

    if (!THB)
    {
        GISMO_ASSERT((dynamic_cast<gsTensorBSpline<2,real_t>*>(&mp.patch(0)) != nullptr),
                     "The input geometry is not a tensor B-spline.");
        gsTensorBSpline<2,real_t> & patch = static_cast<gsTensorBSpline<2,real_t>&>(mp.patch(0));
        patch.degreeElevate(numElev);
        for (index_t i=0; i!=numUHRef; ++i)
            patch.uniformRefine(1);


        for (index_t i=0; i!=patch.domainDim(); ++i)
        {
            index_t size = math::ceil((numHRef+1)*(bbox(i,1) - bbox(i,0)) / l0);
            real_t h=1./size;
            real_t knot = h;
            for (index_t j=0; j<size-1; ++j)
            {
                patch.insertKnot(knot,i);
                knot += h;
            }
        }
        gsInfo<<"New basis = "<<patch.basis()<<"\n";

        gsWrite(mp,outputdir+"result");

        // gsDebugVar(sizes);

    }
    else
    {
        GISMO_ASSERT((dynamic_cast<gsTensorBSpline<2,real_t>*>(&mp.patch(0)) != nullptr),
                     "The input geometry is not a tensor B-spline.");
        gsTensorBSpline<2,real_t> & patch = static_cast<gsTensorBSpline<2,real_t>&>(mp.patch(0));
        gsMatrix<> bbox;
        mp.boundingBox(bbox);
        gsVector<> size = bbox.col(1) - bbox.col(0);
        real_t L = size[0];
        real_t H = size[1];
        if (L > H)
            patch.uniformRefine(math::ceil(L/H)-1,1,0);
        else if (L < H)
            patch.uniformRefine(math::ceil(H/L)-1,1,1);

        for (index_t i=0; i!=numUHRef; ++i)
            patch.uniformRefine(1);

        gsTHBSpline<2,real_t> thb(patch);
        gsMultiPatch<real_t> result;
        result.addPatch(thb);

        gsAdaptiveMeshing<real_t> mesher(result);


        l0 /= H;

        real_t htarget = l0/numHRef;
        real_t hmin = 1;
        for (index_t it=0; it!=10 && hmin>htarget; it++)
        {
            auto domIt  = result.basis(0).domain()->beginAll();
            auto domEnd = result.basis(0).domain()->endAll();
            std::vector<real_t> marked(domEnd-domIt,false);
            hmin = 1;
            for (; domIt<domEnd; ++domIt)
            {
                if (
                    (domIt.centerPoint()[1] >= 0.5-l0 && domIt.centerPoint()[1] <= 0.5+l0) ||
                    (domIt.lowerCorner()[1] >= 0.5-l0 && domIt.lowerCorner()[1] <= 0.5+l0) ||
                    (domIt.upperCorner()[1] >= 0.5-l0 && domIt.upperCorner()[1] <= 0.5+l0)
                   )
                {
                    hmin = math::min(hmin,(domIt.upperCorner()[0] - domIt.lowerCorner()[0])/2);
                    marked[domIt.id()] = true;
                }
            }
            gsInfo<<"Minimum element size: "<<hmin<<", target = "<<htarget<<"\n";

            mesher.options().setSwitch("Admissible",true);
            mesher.options().setInt("RefineRule",1);
            mesher.options().setReal("RefineParam",0.9);
            mesher.options().setInt("MaxLevel",12);
            mesher.getOptions();
            gsHBoxContainer<2,real_t> markedRef;
            mesher.markRef_into(marked,markedRef);
            mesher.refine(markedRef);
            mesher.rebuild();
        }

        gsWrite(result,outputdir+"result");
    }




    return 0;
} // end main

