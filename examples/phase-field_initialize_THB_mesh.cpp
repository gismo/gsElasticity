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

template<short_t dim, class T>
gsMultiPatch<T> createGeometry(const gsMultiPatch<T> & mp)
{
    gsMultiPatch<T> mp_THB;
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
        if ((dynamic_cast<gsTensorBSpline<dim,T>*>(&mp.patch(p))))
        {
            gsTensorBSpline<dim,T> & patch = static_cast<gsTensorBSpline<dim,T>&>(mp.patch(p));
            gsTHBSpline<dim,T> thb(patch);
            mp_THB.addPatch(thb);
        }
        else if ((dynamic_cast<gsTHBSpline<dim,T>*>(&mp.patch(p))))
        {
            gsTHBSpline<dim,T> & patch = static_cast<gsTHBSpline<dim,T>&>(mp.patch(p));
            mp_THB.addPatch(patch);
        }
        else
        {
            GISMO_ERROR("The input geometry is not a tensor B-spline.");
        }
    }
    return mp_THB;
}

template<short_t dim, class T>
void refineGeometry(gsMultiPatch<T> & mp_THB, const gsFunction<T> & crack, gsOptionList mesherOptions)
{
    gsAdaptiveMeshing<dim,T> mesher(mp_THB);
    mesher.options().update(mesherOptions,gsOptionList::ignoreIfUnknown);
    mesher.getOptions();

    gsMatrix<T,dim,2> corners;
    gsMatrix<T> points(dim,math::pow(2,dim)+1);
    points.setZero();
    gsMatrix<T> vals;
    T lowerBound = 0.1;
    T upperBound = 1.0;
    // for (index_t it=0; it!=10 && hmin>htarget; it++)
    for (index_t it=0; it!=mesherOptions.getInt("MaxLevel"); it++)
    {
        gsInfo<<"Refinement iteration "<<it<<":\n";
        gsInfo<<"  Number of elements: "<<mp_THB.basis(0).numElements()<<"\n";
        auto domIt  = mp_THB.basis(0).domain()->beginAll();
        auto domEnd = mp_THB.basis(0).domain()->endAll();
        std::vector<T> marked(domEnd-domIt,false);
        for (; domIt<domEnd; ++domIt)
        {
            // Define the points
            corners.col(0) = domIt.lowerCorner();
            corners.col(1) = domIt.upperCorner();
            gsVector<index_t,dim> np;
            np.setConstant(2);
            gsGridIterator<T,CUBE,dim> grid(corners,np);
            points.col(0) = domIt.centerPoint();
            points.block(0,1,points.rows(),points.cols()-1) = grid.toMatrix();
            crack.piece(0).eval_into(points,vals);
            marked[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
        }
        gsHBoxContainer<dim,T> markedRef;
        mesher.markRef_into(marked,markedRef);
        mesher.refine(markedRef);
        mesher.rebuild();
        gsInfo<<"  Number of elements after refinement: "<<mp_THB.basis(0).numElements()<<"\n";
    }
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numElX = 0;
    index_t numElY = 0;
    index_t numElZ = 0;
    index_t numElev = 0;
    std::string output;
    std::string parInput;
    std::string inputDir;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
    cmd.addInt("x", "numElX","Number of elements in the x direction", numElX);
    cmd.addInt("y", "numElY","Number of elements in the y direction", numElY);
    cmd.addInt("z", "numElZ","Number of elements in the z direction", numElZ);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    inputDir = inputDir + gsFileManager::getNativePathSeparator();
    std::string parInputPath = (parInput.empty() ? inputDir + "parameters.xml" : parInput);
    GISMO_ASSERT(gsFileManager::fileExists(parInputPath), "Input parameter file "<<parInputPath<<" not found.");
    gsInfo << "Input parameter file "<< parInputPath <<"\n";

    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);


    gsFileData<> fd_pars(parInput.empty() ? inputDir + "parameters.xml" : parInput);
    gsInfo << "Input file "<< fd_pars.lastPath() <<"\n";

    GISMO_ASSERT(fd_pars.hasLabel("geometry"), "Geometry not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("crack"), "Crack not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("material"), "Material parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("meshing"), "Adaptive meshing parameters not found in the input file.");

    gsMultiPatch<> mp;
    fd_pars.getLabel("geometry", mp);
    mp.degreeIncrease(numElev);
    mp.uniformRefine(numElX,1,0);
    mp.uniformRefine(numElY,1,1);
    if (mp.geoDim() == 3)
        mp.uniformRefine(numElZ,1,2);

    gsMultiBasis<> mb(mp);
    short_t degree = mb.maxCwiseDegree();

    gsMultiPatch<> crack;
    fd_pars.getLabel("crack", crack);

    gsOptionList materialParameters;
    fd_pars.getLabel("material", materialParameters);

    real_t l0 = materialParameters.getReal("l0");
    real_t Gc = materialParameters.getReal("Gc");
    real_t beta = materialParameters.getReal("beta");

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);

    // gsRBFCurve<real_t, Hat> RBFCurve(crack, beta, beta);
    // Take a range of 2*beta for the mesh!
    gsRBFCurve<real_t, Constant> RBFCurve(crack, (degree)*beta, (degree)*beta);

    gsWriteParaview(mp,RBFCurve,outputdir+"initial",10000);
    gsWriteParaview(mp,outputdir+"mp",10);
    gsWriteParaview(crack,outputdir+"crack",10);

    gsMultiPatch<> mp_THB;
    switch (mp.domainDim())
    {
        case 2:
            mp_THB = createGeometry<2,real_t>(mp);
            refineGeometry<2,real_t>(mp_THB, RBFCurve, mesherOptions);
            break;
        case 3:
            mp_THB = createGeometry<3,real_t>(mp);
            refineGeometry<3,real_t>(mp_THB, RBFCurve, mesherOptions);
            break;
        default:
            GISMO_ERROR("Invalid geometry dimension.");
    }

    gsFileData<> fd_out;
    fd_out.addWithLabel(mp_THB,"geometry");
    fd_out.save(outputdir+"geometry");

    if (plot)
    {
        gsMesh<> mesh(mp_THB.basis(0));
        gsWriteParaview(mesh,outputdir+"THB_mesh",false);
    }



    return 0;
} // end main

