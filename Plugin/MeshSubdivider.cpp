// MeshReducer.C
// Copyright (c) 2009 The Foundry Visionmongers Ltd.  All Rights Reserved.

#include "Subdivider.h"

#include "DDImage/GeoOp.h"
#include "DDImage/SourceGeo.h"
#include "DDImage/ModifyGeo.h"
#include "DDImage/Scene.h"
#include "DDImage/Knob.h"
#include "DDImage/Knobs.h"
#include "DDImage/Mesh.h"
#include "DDImage/Channel3D.h"

#include "DDImage/PolyMesh.h"
#include "DDImage/Triangle.h"

#include <assert.h>


using namespace DD::Image;

static const char* const CLASS = "MeshSubdivider";
static const char* const HELP = "smooth geo faces. triangles and quads";

static const char* const scheme_type_mode[] = { "Bilinear", "Catmull Clark", "Loop", 0 };



class MeshReducer : public SourceGeo
{
private:
    int _enum_scheme;

    Sdc::SchemeType schemeType;
    bool tessQuadsFlag;
    int tessUniformRate;

protected:
    void _validate(bool for_real) override
    {
        // Validate the inputs:
        input0()->validate(for_real);

        // Calculate the geometry hashes:
        GeoOp::_validate(for_real);

        //update_geometry_hashes(); // calls get_geometry_hash()
        set_rebuild(Mask_Geometry);
    }

public:
    static const Description description;
    const char* Class() const { return CLASS; }
    const char* node_help() const { return HELP; }

    MeshReducer(Node* node) : SourceGeo(node)
        , _enum_scheme(1)
        , schemeType(Sdc::SCHEME_CATMARK)
        , tessQuadsFlag(false)
        , tessUniformRate(5)
    {}

    int minimum_inputs() const { return 1; }
    int maximum_inputs() const { return 1; }

    bool test_input(int input, Op* op) const
    {
        //if (input == 1)
        //    return dynamic_cast<AxisOp*>(op) != 0;
        return GeoOp::test_input(input, op);
    }

    Op* default_input(int input) const
    {
        //if (input == 1)
        //    return 0;
        return GeoOp::default_input(input);
    }

    const char* input_label(int input, char* buffer) const
    {
        switch (input) {
        case 0:
            return 0;
        //case 1:
        //    return "axis/cam";
        default:
            return 0;
        }
    }

    void knobs(Knob_Callback f)
    {
        //GeoOp::knobs(f);

        Enumeration_knob(f, &_enum_scheme, scheme_type_mode, "scheme_type");

        Bool_knob(f, &tessQuadsFlag, "preserve_quads"); 
        SetFlags(f, Knob::HIDDEN);

        Int_knob(f, &tessUniformRate, "resolution");

    }

    int knob_changed(Knob* k)
    {
        if (k->is("scheme_type")) {
            switch (_enum_scheme) {
            case 0:

                schemeType = (Sdc::SCHEME_BILINEAR);
                break;
            case 1:

                schemeType = (Sdc::SCHEME_CATMARK);
                break;
            case 2:

                schemeType = (Sdc::SCHEME_LOOP);
                break;
            
            }

            return 1;
        }
        return 1;
    }

    /*! Hash up knobs that affect the primitive attributes. */
    void get_geometry_hash();


    virtual void create_geometry(Scene& scene, GeometryList& out) override;

};

static Op* build(Node* node) { return new MeshReducer(node); }
const Op::Description MeshReducer::description(CLASS, build);

// end of MeshReducer.C

void MeshReducer::get_geometry_hash()
{
    // Get all hashes up-to-date
    GeoOp::get_geometry_hash(); //needs to be GeoOp as it updates the input geo hash
    //SourceGeo::get_geometry_hash(); // probably not necessary

    Hash knob_hash;
    knob_hash.reset();
    knob_hash.append(_enum_scheme);
    knob_hash.append(tessQuadsFlag);
    knob_hash.append(tessUniformRate);
    knob_hash.append(Op::hash()); //not sure..


    // Change the point or vertex attributes hash:
    geo_hash[Group_Primitives].append(knob_hash);
    geo_hash[Group_Points].append(knob_hash);
    geo_hash[Group_Attributes].append(knob_hash);
    geo_hash[Group_Object].append(knob_hash);

}

void MeshReducer::create_geometry(Scene& scene, GeometryList& out)
{
    GeometryList object_list;

    input0()->get_geometry(scene, object_list);

    scene.clear();

    out.delete_objects();

    // Call the engine on all the caches:
    for (uint32_t obj = 0; obj < object_list.objects(); obj++) {
        GeoInfo& info = object_list[obj];

        std::vector<float>     meshVtxPositions;
        std::vector<float>     meshFVarUVs;

        Far::TopologyRefiner* meshTopology = 0;
        meshTopology = readTopologyRefiner(info, this, schemeType, meshVtxPositions, meshFVarUVs);

        //if (meshTopology == NULL) this->error("missing geo"); return;
        //
        
        //  Expand the loaded position and UV arrays to include additional
        //  data (initialized with -1 for distinction), e.g. add a 4-tuple
        //  for RGBA color to the vertex data and add a third field ("w")
        //  to the face-varying data:
        
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
        std::vector<ShapeData> shape_list = tessellateToObj(*meshTopology, vtxData, vtxSize, fvarData, fvarSize, tessQuadsFlag, tessUniformRate);

        const int tessFacetSize = 3 + tessQuadsFlag;

        int vertex_offset = 0;

        for (const auto& shape : shape_list) {
        //for (int i = 50; i < std::min(shape_list.size(), (size_t)150); ++i){
        //    const ShapeData& shape = shape_list[i];
            const uint32_t faces_count = shape.numFacets;
            const size_t vertex_count = shape.outPos.size() / tessFacetSize;

            if (rebuild(Mask_Primitives)) {
                out.add_object(obj);

                auto polymesh = new PolyMesh(faces_count, tessFacetSize);
                auto iterFacets = shape.outFacets.begin();
                for (size_t i = 0; i < faces_count; ++i) {
                    std::vector<int> face; // int face[3] = { *iterFacets++ - vertex_count * obj, *iterFacets++ - vertex_count * obj, *iterFacets++ - vertex_count * obj };
                    face.reserve(tessFacetSize);
                    for (int i = 0; i < tessFacetSize; ++i) face.push_back(*iterFacets++ - vertex_offset);
                    polymesh->add_face(tessFacetSize, face.data());
                }
                out.add_primitive(obj, polymesh);

                // adding material
                out[obj].material = info.material;
            }

            if (rebuild(Mask_Points)) {
                PointList& points = *out.writable_points(obj);
                points.resize(vertex_count);

                auto iterFloat = shape.outPos.begin();

                for (int i = 0; i < vertex_count; ++i) {
                    points[i] = { *iterFloat++, *iterFloat++, *iterFloat++ };
                }
            }


            // Force points and attributes to update:
            set_rebuild(Mask_Points | Mask_Attributes);
            
            if (rebuild(Mask_Attributes)) {
                // Add the UV and N mapping to allow rendering.
                Attribute* uv = out.writable_attribute(obj, Group_Points, "uv", VECTOR4_ATTRIB);
                assert(uv != nullptr);
                uv->resize(shape.outUV.size() / 2);
                unsigned v = 0; //vertex_offset;
                //size_t iter = 0;
                //for (uint32_t i = 0; i < faces_count; ++i) {

                //    for (uint32_t p = 0; p < tessFacetSize; ++p) {
                //        uv->vector4(v++).set(shape.outCoords[outCoordsIter++], shape.outCoords[outCoordsIter++], 0, 1);
                //    }
                //}
                for (uint32_t i = 0; i < vertex_count; ++i)
                    uv->vector4(v++).set(shape.outUV[i * 2], shape.outUV[i * 2 + 1], 0, 1);
            
                //Attribute* N = out.writable_attribute(obj, Group_Vertices, "N", NORMAL_ATTRIB);
                //assert(N != nullptr);
                //N->resize(tri_count * 3);
                //
                //iter = 0;
                //for (uint32_t i = 0; i < tri_count; ++i) {
                //    const auto& normal = mesh.triangles[i].normal;
                //    for (uint32_t p = 0; p < 3; ++p) {
                //        N->normal(iter++).set(normal.x, normal.y, normal.z);
                //    }
                //}
            }

            out.synchronize_objects();

            //info.print_info(std::cout);
            obj++;

            vertex_offset += vertex_count;
        }

        delete meshTopology;
    }

}

