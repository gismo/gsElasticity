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
    std::string output;
    std::string input;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
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

    gsMultiPatch<> mp;
    fd.getLabel("geometry", mp);

    mp.degreeIncrease(numElev);
    mp.uniformRefine(numElX,1,0);
    mp.uniformRefine(numElY,1,1);

    for (size_t p=0; p!=mp.nPatches(); ++p)
        gsInfo<<"Patch "<<p<<": "<<mp.basis(p)<<"\n";

    gsFileData<> fd_out;
    fd_out.addWithLabel(mp,outputdir+"geometry");
    fd_out.save(outputdir+"geometry");

    if (plot) gsWriteParaview(mp,outputdir+"mp",10,true);

    return 0;
} // end main

