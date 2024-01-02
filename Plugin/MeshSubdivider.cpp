// MeshReducer.C
// Copyright (c) 2009 The Foundry Visionmongers Ltd.  All Rights Reserved.

#include "Subdivider.h"

#include "DDImage/GeoOp.h"
#include "DDImage/SourceGeo.h"
#include "DDImage/ModifyGeo.h"
#include "DDImage/Scene.h"
#include "DDImage/Knob.h"
#include "DDImage/Knobs.h"

#include "DDImage/PolyMesh.h"

#include "DDImage/Material.h"

#include <assert.h>


using namespace DD::Image;

static const char* const CLASS = "MeshSubdivider";
static const char* const HELP = "smooth geo faces. into quads only";

static const char* const scheme_type_mode[] = { "Bilinear", "Catmull Clark", /*"Loop",*/ 0 };
static const char* const boundary_interpolation_mode[] = { "boundary_none", "boundary edge only", "boundary edge and corner", 0};

void
getNormal(float N[3], float const du[3], float const dv[3]) {

    N[0] = du[1] * dv[2] - du[2] * dv[1];
    N[1] = du[2] * dv[0] - du[0] * dv[2];
    N[2] = du[0] * dv[1] - du[1] * dv[0];

    float lenSqrd = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
    if (lenSqrd <= 0.0f) {
        N[0] = 0.0f;
        N[1] = 0.0f;
        N[2] = 0.0f;
    }
    else {
        float lenInv = 1.0f / std::sqrt(lenSqrd);
        N[0] *= lenInv;
        N[1] *= lenInv;
        N[2] *= lenInv;
    }
}


class MeshSubdivider : public ModifyGeo
{
private:
    int _enum_scheme;
    int _enum_boundary;

    Sdc::SchemeType schemeType;
    Sdc::Options options;
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

    MeshSubdivider(Node* node) : ModifyGeo(node)
        , _enum_scheme(1)
        , schemeType(Sdc::SCHEME_CATMARK)
        , tessQuadsFlag(false)
        , tessUniformRate(5)
    {}

    int minimum_inputs() const { return 1; }
    int maximum_inputs() const { return 1; }

    //bool test_input(int input, Op* op) const
    //{
    //    //if (input == 0)
    //    //    return dynamic_cast<GeoOp*>(op) != 0;
    //    return GeoOp::test_input(input, op);
    //}

    //Op* default_input(int input) const
    //{
    //    //if (input == 1)
    //    //    return 0;
    //    return GeoOp::default_input(input);
    //}

    //const char* input_label(int input, char* buffer) const
    //{
    //    switch (input) {
    //    //case 0:
    //    //    return "material";
    //    //case 1:
    //    //    return "geo";
    //    default:
    //        return 0;
    //    }
    //}

    void knobs(Knob_Callback f)
    {
        //GeoOp::knobs(f);

        Enumeration_knob(f, &_enum_scheme, scheme_type_mode, "scheme_type");
        Enumeration_knob(f, &_enum_boundary, boundary_interpolation_mode, "boundary_type");

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
        if (k->is("boundary_type")) {
            switch (_enum_boundary) {
            case 0:
                options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_NONE);
                break;
            case 1:
                options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
                break;
            case 2:
                options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
                break;
            }
            return 1;
        }
        return 1;
    }

    /*! Hash up knobs that affect the primitive attributes. */
    void get_geometry_hash();


    virtual void modify_geometry(int obj, Scene& scene, GeometryList& out) override;

};

static Op* build(Node* node) { return new MeshSubdivider(node); }
const Op::Description MeshSubdivider::description(CLASS, build);

// end of MeshReducer.C

void MeshSubdivider::get_geometry_hash()
{
    // Get all hashes up-to-date
    GeoOp::get_geometry_hash(); //needs to be GeoOp as it updates the input geo hash

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

void MeshSubdivider::modify_geometry(int obj, Scene& scene, GeometryList& out)
{
    if (rebuild(Mask_Primitives)) {

        //GeometryList object_list;

        //input0()->get_geometry(scene, out);

        //scene.clear();

        //out.delete_objects();

        // Call the engine on all the caches:
        //for (uint32_t obj = 0; obj < out.objects(); obj++) {
        //    const DD::Image::GeoInfo& info = out[obj];


        //    far_subdivision_with_primvar(info, this, out, obj, std::max(tessUniformRate, 1), schemeType, options);

        //}
        const DD::Image::GeoInfo& info = out[obj];
        far_subdivision_with_primvar(info, this, out, obj, std::max(tessUniformRate, 1), schemeType, options);


        out.synchronize_objects();
        // Force points and attributes to update:
        //set_rebuild(Mask_Points | Mask_Attributes);
    }

}

