/// This is a simple numerical example of modeling active muscle behavior. It is based on the following paper:
/// M.H.Gfrerer and B.Simeon "Fiber-based modeling and simulation of skeletal muscles" 2020
/// Muscle behavior is modeled with the incompressible nonlinear elasticity equations (pressure-displacement formulation)
/// with tendon and muscle materials for the passive part and a fiber-based active response part.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsGeoUtils.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is a simulation of active muscle behavior.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/bicepsRight.xml";

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);

    gsMultiPatch<> muscleMP;
    genMuscleMP(geometry.patch(0),muscleMP);

    return 0;
}
