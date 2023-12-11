//------------------------------------------------------------------------------
//  Tutorial description:
//
//      This tutorial builds on the previous tutorial that makes use of the
//      SurfaceFactory, Surface and Tessellation classes for evaluating and
//      tessellating the limit surface of faces of a mesh by illustrating
//      how the presence of additional data in the mesh arrays is handled.
//
//      As in the previous tutorial, vertex positions and face-varying UVs
//      are provided with the mesh to be evaluated. But here an additional
//      color is interleaved with the position in the vertex data of the
//      mesh and a third component is added to face-varying UV data (making
//      it (u,v,w)).
//
//      To evaluate the position and 2D UVs while avoiding the color and
//      unused third UV coordinate, the Surface::PointDescriptor class is
//      used to describe the size and stride of the desired data to be
//      evaluated in the arrays of mesh data.
//

#include "DDImage/GeoOp.h"
#include "DDImage/GeoInfo.h"

#include <opensubdiv/far/topologyRefiner.h>
#include <opensubdiv/bfr/refinerSurfaceFactory.h>
#include <opensubdiv/bfr/surface.h>
#include <opensubdiv/bfr/tessellation.h>

#include <vector>
#include <string>
#include <cstring>
#include <cstdio>

//  Local headers with support for this tutorial in "namespace tutorial"
//#include "./meshLoader.h"

#include "../../../regression/common/far_utils.h"

//#include "./objWriter.h"



using namespace OpenSubdiv;

struct ShapeData {
    std::vector<float>  outCoords;
    std::vector<float>  outPos, outDu, outDv;
    std::vector<float>  outUV;
    std::vector<int>    outFacets;
    int                 numFacets;
};

//
//  Simple command line arguments to provide input and run-time options:
//
class Args {
public:
    std::string     inputObjFile;
    std::string     outputObjFile;
    Sdc::SchemeType schemeType;
    int             tessUniformRate;
    bool            tessQuadsFlag;
    bool            uv2xyzFlag;

public:
    Args(int argc, char* argv[]) :
        inputObjFile(),
        outputObjFile(),
        schemeType(Sdc::SCHEME_CATMARK),
        tessUniformRate(5),
        tessQuadsFlag(false),
        uv2xyzFlag(false) {

        for (int i = 1; i < argc; ++i) {
            if (strstr(argv[i], ".obj")) {
                if (inputObjFile.empty()) {
                    inputObjFile = std::string(argv[i]);
                }
                else {
                    fprintf(stderr,
                        "Warning: Extra Obj file '%s' ignored\n", argv[i]);
                }
            }
            else if (!strcmp(argv[i], "-o")) {
                if (++i < argc) outputObjFile = std::string(argv[i]);
            }
            else if (!strcmp(argv[i], "-bilinear")) {
                schemeType = Sdc::SCHEME_BILINEAR;
            }
            else if (!strcmp(argv[i], "-catmark")) {
                schemeType = Sdc::SCHEME_CATMARK;
            }
            else if (!strcmp(argv[i], "-loop")) {
                schemeType = Sdc::SCHEME_LOOP;
            }
            else if (!strcmp(argv[i], "-res")) {
                if (++i < argc) tessUniformRate = atoi(argv[i]);
            }
            else if (!strcmp(argv[i], "-quads")) {
                tessQuadsFlag = true;
            }
            else if (!strcmp(argv[i], "-uv2xyz")) {
                uv2xyzFlag = true;
            }
            else {
                fprintf(stderr,
                    "Warning: Unrecognized argument '%s' ignored\n", argv[i]);
            }
        }
    }

private:
    Args() { }
};

Shape* read_GeoOp_input(const DD::Image::GeoInfo& geoInfo, DD::Image::Op* op, Scheme shapescheme = kCatmark, bool isLeftHanded = false);

//
//  The main tessellation function:  given a mesh and vertex positions,
//  tessellate each face -- writing results in Obj format.
//
std::vector<ShapeData>
tessellateToObj(Far::TopologyRefiner const& meshTopology,
    std::vector<float>   const& meshVtxData, int vtxDataSize,
    std::vector<float>   const& meshFVarData, int fvarDataSize,
    const bool                  tessQuadsFlag,
    const int                   tessUniformRate);

//
//  Create a TopologyRefiner from a specified Obj file:
//
Far::TopologyRefiner*
readTopologyRefiner(
    const DD::Image::GeoInfo& geoInfo,
    DD::Image::Op* op,
    Sdc::SchemeType schemeType,
    std::vector<float>& posVector,
    std::vector<float>& uvVector);

////generate GeoList output
//void output_geo(const std::vector<float> outCoords,
//    const std::vector<float> outPos, const std::vector<float> outDu, const std::vector<float> outDv,
//    const std::vector<float> outUV,
//    const std::vector<int>   outFacets,
//    DD::Image::GeometryList& out);

#if 0

//
//  Load command line arguments, specified or default geometry and process:
//
int
main(int argc, char* argv[]) {

    Args args(argc, argv);

    Far::TopologyRefiner* meshTopology = 0;
    std::vector<float>     meshVtxPositions;
    std::vector<float>     meshFVarUVs;


    meshTopology = readTopologyRefiner(
        args.inputObjFile, args.schemeType, meshVtxPositions, meshFVarUVs);
    if (meshTopology == 0) {
        return EXIT_FAILURE;
    }

    //
    //  Expand the loaded position and UV arrays to include additional
    //  data (initialized with -1 for distinction), e.g. add a 4-tuple
    //  for RGBA color to the vertex data and add a third field ("w")
    //  to the face-varying data:
    //
    int numPos = (int)meshVtxPositions.size() / 3;
    int vtxSize = 7;
    std::vector<float> vtxData(numPos * vtxSize, -1.0f);
    for (int i = 0; i < numPos; ++i) {
        vtxData[i * vtxSize] = meshVtxPositions[i * 3];
        vtxData[i * vtxSize + 1] = meshVtxPositions[i * 3 + 1];
        vtxData[i * vtxSize + 2] = meshVtxPositions[i * 3 + 2];
    }

    int numUVs = (int)meshFVarUVs.size() / 2;
    int fvarSize = 3;
    std::vector<float> fvarData(numUVs * fvarSize, -1.0f);
    for (int i = 0; i < numUVs; ++i) {
        fvarData[i * fvarSize] = meshFVarUVs[i * 2];
        fvarData[i * fvarSize + 1] = meshFVarUVs[i * 2 + 1];
    }

    //
    //  Pass the expanded data arrays along with their respective strides:
    //
    tessellateToObj(*meshTopology, vtxData, vtxSize, fvarData, fvarSize, args);

    delete meshTopology;
    return EXIT_SUCCESS;
}

#endif