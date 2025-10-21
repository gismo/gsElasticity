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
    std::vector<index_t> numElX = {0};
    std::vector<index_t> numElY = {0};
    std::vector<index_t> numElZ = {0};
    std::vector<index_t> numElev = {0};
    std::string outputDir;
    std::string parInput;
    std::string inputDir;
    bool into = false;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addMultiInt("e", "numElev","Degree elevation",numElev);
    cmd.addMultiInt("x", "numElX","Number of elements in the x direction", numElX);
    cmd.addMultiInt("y", "numElY","Number of elements in the y direction", numElY);
    cmd.addMultiInt("z", "numElZ","Number of elements in the z direction", numElZ);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "outputDir", "Output directory", outputDir);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);
    cmd.addSwitch("into", "Write the result into the input directory", into);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    char sep = gsFileManager::getNativePathSeparator();
    inputDir = inputDir + gsFileManager::getNativePathSeparator();
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

    gsMultiPatch<> mp;
    fd_pars.getLabel("geometry", mp);
    if (numElev.size()==1)
        numElev = std::vector<index_t>(mp.nPatches(), numElev[0]);
    else
        GISMO_ASSERT(numElev.size()==mp.nPatches(), "Number of entries for degree elevation must be either 1 or equal to the number of patches ("<<mp.nPatches()<<")");
    if (numElX.size()==1)
        numElX  = std::vector<index_t>(mp.nPatches(), numElX[0]);
    else
        GISMO_ASSERT(numElX.size()==mp.nPatches(), "Number of entries for number of elements in X direction must be either 1 or equal to the number of patches ("<<mp.nPatches()<<")");
    if (numElY.size()==1)
        numElY  = std::vector<index_t>(mp.nPatches(), numElY[0]);
    else
        GISMO_ASSERT(numElY.size()==mp.nPatches(), "Number of entries for number of elements in Y direction must be  either 1 or equal to the number of patches ("<<mp.nPatches()<<")");
    if (mp.dim()==3 && numElZ.size()==1)
        numElZ  = std::vector<index_t>(mp.nPatches(), numElZ[0]);
    else if (mp.dim()==3)
        GISMO_ASSERT(numElZ.size()==mp.nPatches(), "Number of entries for number of elements in Y direction must be  either 1 or equal to the number of patches ("<<mp.nPatches()<<")");
    else {}

    GISMO_ENSURE(numElev.size()==mp.nPatches() &&
                 numElX.size()==mp.nPatches() &&
                 numElY.size()==mp.nPatches() &&
                 (mp.dim()!=3 || numElZ.size()==mp.nPatches()), "Arguments -x -y -z must have either one or "<<mp.nPatches()<<" entries.");

    for (size_t p=0; p!=mp.nPatches(); ++p)
    {
        mp.patch(p).degreeIncrease(numElev[p]);
        mp.patch(p).uniformRefine(numElX[p],1,0);
        mp.patch(p).uniformRefine(numElY[p],1,1);
        if (mp.geoDim() == 3)
            mp.patch(p).uniformRefine(numElZ[p],1,2);
    }

    for (size_t p=0; p!=mp.nPatches(); ++p)
        gsInfo<<"Patch "<<p<<": "<<mp.basis(p)<<"\n";

    gsFileData<> fd_out;
    fd_out.addWithLabel(mp,"geometry");
    fd_out.save(outputDir+"geometry");

    if (plot)
    {
        gsMesh<> mesh(mp.basis(0));
        gsWriteParaview(mesh,outputDir+"THB_mesh",false);
    }


    return 0;
} // end main

