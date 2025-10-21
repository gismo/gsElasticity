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
#include <gsHSplines/gsHElementMarker.h>
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
    typedef typename gsHElementHelper<dim,T>::HElementContainer HElementContainer;

    gsHElementMarker<dim,T> marker(mp_THB.basis(0));
    marker.options().update(mesherOptions,gsOptionList::ignoreIfUnknown);

    // gsAdaptiveMeshing<dim,T> mesher(mp_THB);
    // mesher.options().update(mesherOptions,gsOptionList::ignoreIfUnknown);
    // mesher.getOptions();

    gsMatrix<T,dim,2> corners;
    T lowerBound = 0.1;
    T upperBound = 1.0;
    // for (index_t it=0; it!=10 && hmin>htarget; it++)
    gsParaviewCollection refined("markedRef");

    for (index_t it=0; it!=mesherOptions.getInt("MaxLevel"); it++)
    {
        gsInfo<<"Refinement iteration "<<it<<":\n";
        index_t numEl = mp_THB.patch(0).basis().numElements();
        gsInfo<<"  Number of elements: "<<numEl<<"\n";
        // auto domIt  = mp_THB.basis(0).domain()->beginAll();
        // auto domEnd = mp_THB.basis(0).domain()->endAll();
        std::vector<T> marked(numEl,false);
        // for (; domIt<domEnd; ++domIt)
        gsStopwatch timer;
        for (auto & domIt : mp_THB.patch(0).basis().domain()->allElements())
        {
            gsMatrix<T> vals;
            gsMatrix<T> points(dim,math::pow(2,dim)+1);
            // gsMatrix<T> points(dim,1);
            points.setZero();

            // Define the points
            corners.col(0) = domIt.lowerCorner();
            corners.col(1) = domIt.upperCorner();

            gsVector<index_t,dim> np;
            np.setConstant(2);
            gsGridIterator<T,CUBE,dim> grid(corners,np);
            points.col(0) = domIt.centerPoint();
            points.block(0,1,points.rows(),points.cols()-1) = grid.toMatrix();
            mp_THB.piece(0).eval_into(points,vals);
            std::swap(points,vals);
            crack.piece(0).eval_into(points,vals);
            marked[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
        }
        gsInfo<<"  Computing values took "<<timer.stop()<<" seconds.\n";
        // gsHBoxContainer<dim,T> markedRef;
        // mesher.markRef_into(marked,markedRef);
        timer.restart();
        marker.setErrors(marked);
        timer.stop();
        gsInfo<<"  Setting errors took "<<timer.stop()<<" seconds.\n";
        timer.restart();
        HElementContainer markedRef = marker.markRef();
        gsInfo<<"  Marking took "<<timer.stop()<<" seconds.\n";
        timer.restart();
        std::vector<index_t> refBox = marker.toRefBoxes(markedRef);
        gsInfo<<"  Conversion to refinement boxes took "<<timer.stop()<<" seconds.\n";
        timer.restart();
        mp_THB.patch(0).refineElements(refBox);
        gsInfo<<"  Refinement took "<<timer.stop()<<" seconds.\n";
        gsInfo<<"  Number of elements after refinement: "<<mp_THB.basis(0).numElements()<<"\n";


        ///////////////////////////////////////////////////////////////////////////////////////////
        // PLOT
        ///////////////////////////////////////////////////////////////////////////////////////////
        gsMatrix<T> boxes;
        gsVector<size_t> levels;
        std::tie(boxes,levels) = marker.helper().toBoxesAndLevels(markedRef);
        gsWriteParaview(boxes,"markedRef_"+util::to_string(it),gsVector<real_t>(levels.cast<real_t>()));
        refined.addPart("markedRef_"+util::to_string(it)+".vtu",it,"Solution");
    }
    refined.save();
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool plotMesh = false;
    index_t numElX = 0;
    index_t numElY = 0;
    index_t numElZ = 0;
    index_t numElev = 0;
    std::string outputDir;
    std::string parInput;
    std::string inputDir;
    bool into = false;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
    cmd.addInt("x", "numElX","Number of elements in the x direction", numElX);
    cmd.addInt("y", "numElY","Number of elements in the y direction", numElY);
    cmd.addInt("z", "numElZ","Number of elements in the z direction", numElZ);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("plotMesh","Create a ParaView visualization file with the mesh", plotMesh);
    cmd.addString("o", "outputDir", "Output directory", outputDir);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);
    cmd.addSwitch("into", "Write the result into the input directory", into);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    char sep = gsFileManager::getNativePathSeparator();
    inputDir = inputDir + sep;
    std::string parInputPath = (parInput.empty() ? inputDir + "parameters.xml" : parInput);
    GISMO_ASSERT(gsFileManager::fileExists(parInputPath), "Input parameter file "<<parInputPath<<" not found.");
    gsInfo << "Input parameter file "<< parInputPath <<"\n";

    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    if (into)
        outputDir = inputDir;
    else
    {
        if (outputDir.empty())
            outputDir = std::string(".") + sep + "output" + sep;
        else
            outputDir += sep;
        gsFileManager::mkdir(outputDir);
    }
    gsInfo<< "Output directory: "<<outputDir<<"\n";

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

    gsMatrix<> supp;
    if (fd_pars.hasLabel("support"))
    {
        fd_pars.getLabel("support", supp);
        GISMO_ASSERT(supp.cols() == 2 && supp.rows() == mp.geoDim(),
                        "Support must be a matrix of size "<<mp.geoDim()<<"x2.");
    }
    else
        supp = mp.patch(0).support();

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

    gsMultiPatch<> mp_THB;
    real_t hmin;
    switch (mp.domainDim())
    {
        case 2:
            mp_THB = createGeometry<2,real_t>(mp);
            hmin = (static_cast<gsTHBSplineBasis<2,real_t>&>(mp_THB.basis(0)).tensorLevel(mesherOptions.getInt("MaxLevel")).getMinCellLength());
            break;
        case 3:
            mp_THB = createGeometry<3,real_t>(mp);
            hmin = (static_cast<gsTHBSplineBasis<2,real_t>&>(mp_THB.basis(0)).tensorLevel(mesherOptions.getInt("MaxLevel")).getMinCellLength());
            break;
        default:
            GISMO_ERROR("Invalid geometry dimension.");
    }

    gsInfo<<"Beta = "<<beta<<", hmin = "<<hmin<<", (p+1)*hmin = "<<(degree+1)*hmin<<"\n";
    beta = math::max(beta,hmin);
    beta*= degree+1;
    gsInfo<<"Final beta = "<<beta<<"\n";



    // gsRBFCurve<real_t, Hat> RBFCurve(crack, beta, beta);
    // Take a range of 2*beta for the mesh!
    gsRBFCurve<real_t, Constant> RBFCurve(crack, beta, beta);

    if(plot)
    {
        // Evaluate the geometry in the support
        gsVector<unsigned> npts = uniformSampleCount<real_t>(supp.col(0), supp.col(1), 1000000);
        gsMatrix<> points = gsPointGrid<real_t>(supp.col(0),supp.col(1),npts);
        gsMatrix<> eval_geo, eval_damage;
        mp.piece(0).eval_into(points, eval_geo);
        RBFCurve.piece(0).eval_into(eval_geo, eval_damage);
        gsWriteParaviewTPgrid(eval_geo,eval_damage,npts.template cast<index_t>(),outputDir+"initial");
    }

    gsWriteParaview(mp,outputDir+"mp",10);
    gsWriteParaview(crack,outputDir+"crack",10);

    switch (mp.domainDim())
    {
        case 2:
            refineGeometry<2,real_t>(mp_THB, RBFCurve, mesherOptions);
            break;
        case 3:
            refineGeometry<3,real_t>(mp_THB, RBFCurve, mesherOptions);
            break;
        default:
            GISMO_ERROR("Invalid geometry dimension.");
    }

    gsFileData<> fd_out;
    fd_out.addWithLabel(mp_THB,"geometry");
    fd_out.save(outputDir+"geometry");

    if (plotMesh)
    {
        gsMesh<> mesh(mp_THB.basis(0));
        mp_THB.patch(0).evaluateMesh(mesh);
        gsWriteParaview(mesh,outputDir+"THB_mesh",false);
    }



    return 0;
} // end main

