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
#include <gsModeling/gsRBFCurve.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numElX = 0;
    index_t numElY = 0;
    index_t numElev = 0;
    index_t numHRef = 0;
    std::string output;
    std::string input;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
    cmd.addInt("r", "numHRef","Number of elements in the crack size", numHRef);
    cmd.addInt("x", "numElX","Number of elements in the x direction", numElX);
    cmd.addInt("y", "numElY","Number of elements in the y direction", numElY);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "input", "Input XML file", input);
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

    GISMO_ASSERT(fd.hasLabel("geometry"), "Geometry not found in the input file.");
    GISMO_ASSERT(fd.hasLabel("crack"), "Crack not found in the input file.");
    GISMO_ASSERT(fd.hasLabel("material"), "Material parameters not found in the input file.");

    gsMultiPatch<> mp;
    fd.getLabel("geometry", mp);
    mp.degreeIncrease(numElev);
    mp.uniformRefine(numElX,1,0);
    mp.uniformRefine(numElY,1,1);

    gsMultiPatch<> crack;
    fd.getLabel("crack", crack);

    gsOptionList materialParameters;
    fd.getLabel("material", materialParameters);

    real_t l0 = materialParameters.getReal("l0");
    real_t Gc = materialParameters.getReal("Gc");
    real_t beta = materialParameters.getReal("beta");

    gsRBFCurve<real_t, Constant> RBFCurve(crack, beta, beta);


    gsWriteParaview(mp,RBFCurve,outputdir+"initial",10000);
    gsWriteParaview(mp,outputdir+"mp",10);
    gsWriteParaview(crack,outputdir+"crack",10);

    // Handle the mesh
    gsMultiPatch<> mp_THB;
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
        if ((dynamic_cast<gsTensorBSpline<2,real_t>*>(&mp.patch(p))))
        {
            gsTensorBSpline<2,real_t> & patch = static_cast<gsTensorBSpline<2,real_t>&>(mp.patch(p));
            gsTHBSpline<2,real_t> thb(patch);
            mp_THB.addPatch(thb);
        }
        else if ((dynamic_cast<gsTHBSpline<2,real_t>*>(&mp.patch(p))))
        {
            gsTHBSpline<2,real_t> & patch = static_cast<gsTHBSpline<2,real_t>&>(mp.patch(p));
            mp_THB.addPatch(patch);
        }
        else
        {
            GISMO_ERROR("The input geometry is not a tensor B-spline.");
        }
    }

    gsAdaptiveMeshing<2,real_t> mesher(mp_THB);

    real_t htarget = l0/numHRef;
    real_t hmin = 1;
    gsMatrix<real_t,2,2> corners;
    gsMatrix<real_t> points(2,5);
    gsMatrix<real_t> vals;
    real_t lowerBound = 0.1;
    real_t upperBound = 1.0;
    for (index_t it=0; it!=10 && hmin>htarget; it++)
    {
        auto domIt  = mp_THB.basis(0).domain()->beginAll();
        auto domEnd = mp_THB.basis(0).domain()->endAll();
        std::vector<real_t> marked(domEnd-domIt,false);
        hmin = 1;
        for (; domIt<domEnd; ++domIt)
        {
            // Define the points
            corners.col(0) = domIt.lowerCorner();
            corners.col(1) = domIt.upperCorner();
            gsGridIterator<real_t,CUBE,2> grid(corners,2);
            points.col(0) = domIt.centerPoint();
            points.block(0,1,2,4) = grid.toMatrix();
            RBFCurve.piece(0).eval_into(points,vals);
            // gsDebugVar(points);
            // gsDebugVar(vals);
            marked[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
            hmin = math::min(hmin,(domIt.upperCorner()[0] - domIt.lowerCorner()[0])/2);
        }
        gsInfo<<"Minimum element size: "<<hmin<<", target = "<<htarget<<"\n";
        // if (gsAsVector<real_t>(marked).sum()==0)
        // {
        //     gsInfo<<"No elements marked for refinement. Stopping refinement.\n";
        //     break;
        // }
        // gsDebugVar(gsAsVector<real_t>(marked).sum());
        // gsDebugVar(gsAsVector<real_t>(marked).size());

        // NOTE: If there are no marked elements, a uniform refinement is performed
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

    gsFileData<> fd_out(outputdir+"geometry");
    fd_out.addWithLabel(mp_THB,outputdir+"geometry");

    gsWriteParaview(mp_THB,outputdir+"THB",10000,true);



    return 0;
} // end main

