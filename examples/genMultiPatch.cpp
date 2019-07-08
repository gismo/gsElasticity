#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string filename = "";

    gsCmdLine cmd("");
    cmd.addPlainString("name","File with boundary curves.",filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsMultiPatch<> bdry;
    gsReadFile<>(filename,bdry);
    gsMultiPatch<> geo;

    for (index_t i = 0; i < bdry.nPatches()/4; i++)
    {
        gsMultiPatch<> patch;
        for (index_t j = 0; j < 4; j++)
            patch.addPatch(bdry.patch(4*i+j));
        gsCoonsPatch<real_t> coonsPatch(patch);
        coonsPatch.compute();
        geo.addPatch(coonsPatch.result());
    }

    geo.computeTopology();
    filename = filename.substr(filename.find_last_of("/\\")+1); // file name without a path
    filename = filename.substr(0,filename.find_last_of(".\\"));
    gsWrite(geo,filename + "_2D");
    return 0;
}
