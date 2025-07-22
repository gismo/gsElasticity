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
    mp.degreeIncrease(numElev);
    mp.uniformRefine(numElX,1,0);
    mp.uniformRefine(numElY,1,1);
    if (mp.geoDim() == 3)
        mp.uniformRefine(numElZ,1,2);

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

