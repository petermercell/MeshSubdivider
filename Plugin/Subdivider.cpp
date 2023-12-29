#include "Subdivider.h"

#include "DDImage/GeoOp.h"
#include "DDImage/PolyMesh.h"

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>

#include <opensubdiv/osd/cpuEvaluator.h>
#include <opensubdiv/osd/cpuVertexBuffer.h>


void far_subdivision_with_primvar(const DD::Image::GeoInfo& geoInfo, DD::Image::GeoOp* op, DD::Image::GeometryList& out, const int obj, const int maxlevel, Sdc::SchemeType shapescheme) {


    enum NormalApproximation normalApproximation = CrossTriangle;

    // 'vertex' primitive variable data & topology
    std::vector< float >g_verts; // [8] [3]

    int g_nverts = 0; //g_verts.size();
    int  g_nfaces = 0; //6

    std::vector< int > g_vertsperface; // [6]

    std::vector< int> g_vertIndices; // [24] 

    // 'face-varying' primitive variable data & topology for UVs
    std::vector< float> g_uvs; // [14] [2]

    int g_nuvs = 0; //g_uvs.size();

    std::vector< int> g_uvIndices; // [24]

    // 'face-varying' primitive variable data & topology for color
    //std::vector< DD::Image::Vector4 >g_colors; // [24] [4]

    //int g_ncolors = 0; //g_colors.size();

    //std::vector< int> g_colorIndices; // [24]

    const DD::Image::PointList* geoPoints = geoInfo.point_list();

    // reserve mesh elements
    uint32_t num_vertices = (uint32_t)geoPoints->size();;

#ifdef DEBUG_MESH
    std::cout << "readgeo --> total vertex amount: " << num_vertices << std::endl;
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------> , UVs, Normals(disabled) vertex
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////



    //////get normal information
    //const DD::Image::AttribContext* N_ref = geoInfo.get_attribcontext("N");
    //const DD::Image::AttributePtr vertex_normal = N_ref ? N_ref->attribute : DD::Image::AttributePtr();

    //if (!vertex_normal) {
    //    op->error("Missing \"N\" channel from geometry");
    //    return false;
    //}

    //////////////////////////////////////////////////////////////////////// adding all normals positions 
    //if (vertex_normal) {

    //    uint32_t numNormals = vertex_normal->size();
    //    s->normals.reserve(numNormals * 3);

    //    for (int iterN = 0; iterN < numNormals; iterN++) {

    //        const DD::Image::Vector3& localN = vertex_normal->normal(iterN);
    //        DD::Image::Vector4 worldN = geoInfo.matrix * localN;

    //        //std::cout << "vertex N values at index " << iterN << ": --> " << localN.x << "  " << localN.y << "  " << localN.z << std::endl;

    //        s->normals.push_back(worldN.x);
    //        s->normals.push_back(worldN.y);
    //        s->normals.push_back(worldN.z);

    //    }
    //}
    //DD::Image::GroupType n_group_type = DD::Image::Group_None;
    //if (vertex_normal)
    //    DD::Image::GroupType n_group_type = N_ref->group;    // normal group type 

    // get the original uv attribute used to restore untouched uv coordinate
    const DD::Image::AttribContext* UV_ref = geoInfo.get_attribcontext("uv");
    DD::Image::AttributePtr uv_original = UV_ref ? UV_ref->attribute : DD::Image::AttributePtr();
    DD::Image::GroupType t_group_type = DD::Image::Group_None;

    //if (!uv_original) {
    //    op->error("Missing uv channel from geometry");
    //    //return false;
    //}
    /*else {*/
    if (uv_original) {

        uint32_t numUVs = uv_original->size();

        t_group_type = UV_ref->group;

        //////////////////////////////////////////////////////////////////////// adding all uv positions 
        if (uv_original)
        {
            g_uvs.reserve(numUVs*2);

            for (uint32_t iterUV = 0; iterUV < numUVs; iterUV++) {

                //std::cout << "input UV at iter " << iterUV << ": --> " << uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w << "  " << uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w << std::endl;

                //g_uvs.push_back({ uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w, uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w });
                g_uvs.push_back(uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w);
                g_uvs.push_back(uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w);
            }
        }
    }

    g_nuvs = g_uvs.size();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////////////////////// ADDING ALL VERTEX POSITIONS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    g_verts.reserve(num_vertices*3);

    for (uint32_t iter = 0; iter < num_vertices; iter++)
    {
        const DD::Image::Vector3& local_point = geoInfo.point_array()[iter];

        DD::Image::Vector4 w_point(geoInfo.matrix * local_point);

        //g_verts.push_back({ w_point[0],w_point[1],w_point[2] });
        g_verts.push_back(w_point[0]);
        g_verts.push_back(w_point[1]);
        g_verts.push_back(w_point[2]);

    }

    g_nverts = g_verts.size();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////--------------------------------------------------------------------> GENERATE MESH
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///tri or quads primitives
    const DD::Image::Primitive** prim_array = geoInfo.primitive_array();
    if (prim_array)
    {
        const int num_prims = geoInfo.primitives(); //amount of primitives


        // vertex index array
        std::vector<uint32_t> varray;

        int vertexCount = 0; // faces * vertex ----- for uvs index
        //int count_for_normals = 0; // faces count for normals smoothing

        for (int p = 0; p < num_prims; ++p) {

            // get primitive
            const DD::Image::Primitive* prim = geoInfo.primitive(p);
            const uint32_t num_faces = prim->faces();


            // iterate all faces in the primitive
            for (uint32_t f = 0; f < num_faces; f++) {

                // get the number of vertex used by face f
                const uint32_t num_verts = prim->face_vertices(f);

                g_nfaces++;
                g_vertsperface.push_back(num_verts);

                // get all used vertices index
                varray.resize(num_verts);

                prim->get_face_vertices(f, &varray[0]);

                // check all vertices in the face
                for (unsigned v = 0; v < num_verts; ++v) {

                    // TODO: vertex_offset is always 0 for polymeshes, since the offset is dependent upon the face in this case
                    // (each subface can have a differing number of vertices). For all other primitive type, the offset is fixed 
                    // per primitive.
                    // Really the Primitive class needs extending to take this into account (vertex_offset() could take a face index 
                    // to achieve this). For now fix this here by using a locally calculated offset.
                    mFnAssert((prim->getPrimitiveType() != ePolyMesh) || (prim->vertex_offset() == 0));
                    const unsigned int vertexOffset = (prim->getPrimitiveType() == DD::Image::ePolyMesh) ? vertexCount : prim->vertex_offset();

                    unsigned vi = prim->vertex(varray[v]); // +prim->vertex_offset();
                    unsigned ti = (t_group_type == DD::Image::Group_Points) ? vi : (vertexOffset + v);
                    //unsigned ni = (n_group_type == DD::Image::Group_Points) ? vi : (vertexOffset + v);

                    g_vertIndices.push_back(vi);
                    if (uv_original)g_uvIndices.push_back(ti);
                    //if (vertex_normal) s->facenormals.push_back(ni);

                }

                vertexCount += num_verts;

            }
        }

    }

    // Populate a topology descriptor with our raw data
    Descriptor desc;
    desc.numVertices = g_nverts;
    desc.numFaces = g_nfaces;
    desc.numVertsPerFace = g_vertsperface.data();
    desc.vertIndicesPerFace = g_vertIndices.data();

    // Create a face-varying channel descriptor
    const int numChannels = 1;
    const int channelUV = 0;
    //const int channelColor = 1;
    Descriptor::FVarChannel channels[numChannels];
    channels[channelUV].numValues = g_nuvs;
    channels[channelUV].valueIndices = g_uvIndices.data();
    //channels[channelColor].numValues = g_ncolors;
    //channels[channelColor].valueIndices = g_colorIndices.data();

    // Add the channel topology to the main descriptor
    desc.numFVarChannels = numChannels;
    desc.fvarChannels = channels;

    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_NONE);
    options.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_ALL);
    //options.SetCreasingMethod(Sdc::Options::CREASE_CHAIKIN);
    //options.SetTriangleSubdivision(Sdc::Options::TRI_SUB_SMOOTH);



    // Instantiate a Far::TopologyRefiner from the descriptor
    Far::TopologyRefiner* refinerOSD =
        Far::TopologyRefinerFactory<Descriptor>::Create(desc,
            Far::TopologyRefinerFactory<Descriptor>::Options(shapescheme, options));

    // Uniformly refine the topolgy up to 'maxlevel'
    // note: fullTopologyInLastLevel must be true to work with face-varying data
    {
        Far::TopologyRefiner::UniformOptions refineOptions(maxlevel);
        refineOptions.fullTopologyInLastLevel = true;
        refinerOSD->RefineUniform(refineOptions);
        //refiner->RefineAdaptive(refineOptions);
    }

    int nCoarseVerts = 0;
    int nRefinedVerts = 0;
  
    //
    // Setup phase
    //
    Far::StencilTable const* stencilTableVTX = NULL;
    { 
        Far::StencilTableFactory::Options options;
        options.generateOffsets = true;
        options.generateIntermediateLevels = false;

        options.interpolationMode = Far::StencilTableFactory::INTERPOLATE_VERTEX;

        stencilTableVTX = Far::StencilTableFactory::Create(*refinerOSD, options);

        nCoarseVerts = refinerOSD->GetLevel(0).GetNumVertices();
        nRefinedVerts = stencilTableVTX->GetNumStencils();

    }

    // Setup a buffer for vertex primvar data:
    Osd::CpuVertexBuffer* vbuffer =
    Osd::CpuVertexBuffer::Create(3, nCoarseVerts + nRefinedVerts);

    //
    // Execution phase (every frame)
    //

    {
        // Pack the control vertex data at the start of the vertex buffer
        // and update every time control data changes
        vbuffer->UpdateData(g_verts.data(), 0, nCoarseVerts);


        Osd::BufferDescriptor srcDesc(0, 3, 3);
        Osd::BufferDescriptor dstDesc(nCoarseVerts * 3, 3, 3);

        // Launch the computation
        Osd::CpuEvaluator::EvalStencils(vbuffer, srcDesc,
            vbuffer, dstDesc,
            stencilTableVTX);
    }

    { 

        printf("vertices ");
        float const* refinedVerts = vbuffer->BindCpuBuffer() + 3 * nCoarseVerts;
        for (int i = 0; i < nRefinedVerts; ++i) {
            float const* vert = refinedVerts + 3 * i;
            printf("-p %f %f %f\n", vert[0], vert[1], vert[2]);
        }
        printf("-c 1;\n");
    }

    int nCoarseUVs = 0;
    int nRefinedUVs = 0;

    //
    // Setup phase
    //
    Far::StencilTable const* stencilTableUV = NULL;
    { 
        Far::StencilTableFactory::Options options;
        options.generateOffsets = true;
        options.generateIntermediateLevels = false;

        options.interpolationMode = Far::StencilTableFactory::INTERPOLATE_FACE_VARYING;

        stencilTableUV = Far::StencilTableFactory::Create(*refinerOSD, options);

        nCoarseUVs = refinerOSD->GetLevel(0).GetNumFVarValues(channelUV);
        nRefinedUVs = stencilTableUV->GetNumStencils();

    }

    // Setup a buffer for vertex primvar data:
    Osd::CpuVertexBuffer* uvbuffer =
        Osd::CpuVertexBuffer::Create(2, nCoarseUVs + nRefinedUVs);

    //
    // Execution phase (every frame)
    //

    {
        // Pack the control vertex data at the start of the vertex buffer
        // and update every time control data changes
        uvbuffer->UpdateData(g_uvs.data(), 0, nCoarseUVs);


        Osd::BufferDescriptor srcDesc(0, 2, 2);
        Osd::BufferDescriptor dstDesc(nCoarseUVs * 2, 2, 2);

        // Launch the computation
        Osd::CpuEvaluator::EvalStencils(uvbuffer, srcDesc,
            uvbuffer, dstDesc,
            stencilTableUV);
    }


    std::vector<UV> uvs_array;
    uvs_array.reserve(nRefinedUVs);

    { // Visualization with Maya : print a MEL script that generates particles
    // at the location of the refined vertices

        printf("uvs ");
        float const* refinedUVs = uvbuffer->BindCpuBuffer() + 2 * nCoarseUVs;
        for (int i = 0; i < nRefinedUVs; ++i) {
            float const* uv = refinedUVs + 2 * i;
            printf("-index %d: %f %f\n", i, uv[0], uv[1]);
            uvs_array.push_back({ uv[0], uv[1] });
        }
        printf("-c 1;\n");
    }

    Far::TopologyLevel const& refined_refLastLevel = refinerOSD->GetLevel(maxlevel);
    int refined_nverts = refined_refLastLevel.GetNumVertices();
    int refined_nfaces = refined_refLastLevel.GetNumFaces();

    if (op->rebuild(Mask_Primitives)) {
        out.add_object(obj);
        //out.writable_primitive(obj, );

        auto polymesh = new PolyMesh(refined_nfaces, 4);
        for (size_t f = 0; f < refined_nfaces; ++f) {
            Far::ConstIndexArray fverts = refined_refLastLevel.GetFaceVertices(f);
            std::vector<int> face;
            face.reserve(4);
            for (int i = 0; i < 4; ++i) {
                //printf("prim: %d, index: %d \n", f, fverts[i]);
                face.push_back(fverts[i] /*- vertex_offset*/);
            }
            polymesh->add_face(4, face.data());
        }
        out.add_primitive(obj, polymesh);

        //// adding material
        //auto input_shader = (Iop*)input0();
        out[obj].material = geoInfo.material;
        //out[obj + 2].useMaterialContext = true;

        // Force points and attributes to update:
        op->set_rebuild(Mask_Points | Mask_Attributes);
    }

    if (op->rebuild(Mask_Points)) {
        PointList& points = *out.writable_points(obj);
        points.resize(refined_nverts);

        float const* refinedVerts = vbuffer->BindCpuBuffer() + 3 * nCoarseVerts;
        for (int i = 0; i < nRefinedVerts; ++i) {
            float const* vert = refinedVerts + 3 * i;
            points[i] = { vert[0], vert[1], vert[2] };
        }

    }

    int nuvs = refined_refLastLevel.GetNumFVarValues(channelUV);
    int firstOfLastUvs = refinerOSD->GetNumFVarValuesTotal(channelUV) - nuvs;
    std::cout << firstOfLastUvs << "\n";
    if (op->rebuild(Mask_Attributes) && nuvs != 0) {
        Attribute* uv = out.writable_attribute(obj, Group_Vertices, "uv", VECTOR4_ATTRIB);
        assert(uv != nullptr);
        uv->resize(nuvs * refined_nfaces);

        int faceCount = 0;
        for (int face = 0; face < refined_nfaces; ++face) {
            Far::ConstIndexArray fverts = refined_refLastLevel.GetFaceVertices(face);
            Far::ConstIndexArray fuvs = refined_refLastLevel.GetFaceFVarValues(face, channelUV);

            // all refined Catmark faces should be quads
            assert(fverts.size() == 4 && fuvs.size() == 4);


            float const* refinedUVs = uvbuffer->BindCpuBuffer() + 2 * nCoarseUVs;
            for (int fvvert = 0; fvvert < fverts.size(); ++fvvert) {
                float const* uvVal = uvs_array[fuvs[fvvert]].GetData();
                printf("uv_out at prim: %d, index: %d, uvs:  %f %f\n", face, fuvs[fvvert], uvVal[0], uvVal[1]);
                uv->vector4(faceCount + fvvert).set(uvVal[0], uvVal[1], 0, 1);
#ifdef  DEBUG
                printf("uv_out %d\n", fuvs[fvvert]);
                printf("uv_out %f %f\n", uvIn[0], uvIn[1]);
#endif
            }     

            faceCount += 4;
        }
    }

    //if (op->rebuild(Mask_Attributes)) {
    //    Attribute* N = out.writable_attribute(obj, Group_Points, "N", NORMAL_ATTRIB);
    //    assert(N != nullptr);
    //    N->resize(nverts);


    //    for (int vert = 0; vert < nverts; ++vert) {
    //        float const* pos = normals[vert].GetPosition();
    //        N->normal(vert).set(pos[0], pos[1], pos[2]);
    //    }
    //}

    delete refinerOSD;

}


