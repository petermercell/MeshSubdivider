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

//#include <opensubdiv/far/topologyRefiner.h>
//#include <opensubdiv/bfr/refinerSurfaceFactory.h>
//#include <opensubdiv/bfr/surface.h>
//#include <opensubdiv/bfr/tessellation.h>

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>

#include <vector>
#include <string>
#include <cstring>
#include <cstdio>

using namespace DD::Image;
using namespace OpenSubdiv;

//#define DEBUG

//------------------------------------------------------------------------------
// Math helpers.
//
//

// Returns the normalized version of the input vector
inline void
normalize(float* n) {
    float rn = 1.0f / sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] *= rn;
    n[1] *= rn;
    n[2] *= rn;
}

// Returns the cross product of \p v1 and \p v2.                                
inline void cross(float const* v1, float const* v2, float* vOut)
{
    vOut[0] = v1[1] * v2[2] - v1[2] * v2[1];
    vOut[1] = v1[2] * v2[0] - v1[0] * v2[2];
    vOut[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

//------------------------------------------------------------------------------
// Face-varying implementation.
//
//
struct Vertex {

    // Minimal required interface ----------------------
    Vertex() {
        Clear();
    }

    Vertex(Vertex const& src) {
        position[0] = src.position[0];
        position[1] = src.position[1];
        position[2] = src.position[2];
    }

    void Clear() {
        position[0] = position[1] = position[2] = 0.0f;
    }

    void AddWithWeight(Vertex const& src, float weight) {
        position[0] += weight * src.position[0];
        position[1] += weight * src.position[1];
        position[2] += weight * src.position[2];
    }

    // Public interface ------------------------------------
    void SetPosition(float x, float y, float z) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
    }

    const float* GetPosition() const {
        return position;
    }

    float* GetPosition() {
        return position;
    }

    // Public interface ------------------------------------
    float const* GetData() const {
        return position;
    }

    float position[3];
};

struct UV {

    // Minimal required interface ----------------------
    UV() {
        Clear();
    }

    UV(UV const& src) {
        uv[0] = src.uv[0];
        uv[1] = src.uv[1];
    }

    UV(const float u, const float v) {
        uv[0] = u;
        uv[1] = v;
    }

    void Clear() {
        uv[0] = uv[1] = uv[2] = 0.0f;
    }

    void AddWithWeight(UV const& src, float weight) {
        uv[0] += weight * src.uv[0];
        uv[1] += weight * src.uv[1];
    }

    // Public interface ------------------------------------
    void Setuv(float x, float y, float z) {
        uv[0] = x;
        uv[1] = y;
    }

    const float* GetPosition() const {
        return uv;
    }

    // Public interface ------------------------------------
    float const* GetData() const {
        return uv;
    }

    float uv[2];
};

//------------------------------------------------------------------------------
// Face-varying container implementation.
//
// We are using a uv texture layout as a 'face-varying' primtiive variable
// attribute. Because face-varying data is specified 'per-face-per-vertex',
// we cannot use the same container that we use for 'vertex' or 'varying'
// data. We specify a new container, which only carries (u,v) coordinates.
// Similarly to our 'Vertex' container, we add a minimaliztic interpolation
// interface with a 'Clear()' and 'AddWithWeight()' methods.
//
struct FVarVertexUV {

    // Minimal required interface ----------------------
    void Clear() {
        u = v = 0.0f;
    }

    void AddWithWeight(FVarVertexUV const& src, float weight) {
        u += weight * src.u;
        v += weight * src.v;
    }

    // Basic 'uv' layout channel
    float u, v;
};

struct FVarVertexColor {

    // Minimal required interface ----------------------
    void Clear() {
        r = g = b = a = 0.0f;
    }

    void AddWithWeight(FVarVertexColor const& src, float weight) {
        r += weight * src.r;
        g += weight * src.g;
        b += weight * src.b;
        a += weight * src.a;
    }

    // Basic 'color' layout channel
    float r, g, b, a;
};

// Approximation methods for smooth normal computations
enum NormalApproximation
{
    CrossTriangle,
    CrossQuad,
    Limit
};

typedef Far::TopologyDescriptor Descriptor;

void far_subdivision_with_primvar(  const DD::Image::GeoInfo& geoInfo,
                                    DD::Image::GeoOp* op,
                                    DD::Image::GeometryList& out,
                                    const int obj,
                                    const int maxlevel,
                                    Sdc::SchemeType shapescheme);
