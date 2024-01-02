#include "Subdivider.h"

#include "DDImage/GeoOp.h"
#include "DDImage/PolyMesh.h"


void far_subdivision_with_primvar(const DD::Image::GeoInfo& geoInfo, DD::Image::GeoOp* op, DD::Image::GeometryList& out, const int obj, const int maxlevel, Sdc::SchemeType shapescheme, Sdc::Options opt) {


    enum NormalApproximation normalApproximation = CrossTriangle;

    // 'vertex' primitive variable data & topology
    std::vector< DD::Image::Vector3 >g_verts; // [8] [3]

    int g_nverts = 0; //g_verts.size();
    int  g_nfaces = 0; //6

    std::vector< int > g_vertsperface; // [6]

    std::vector< int> g_vertIndices; // [24] 

    // 'face-varying' primitive variable data & topology for UVs
    std::vector< DD::Image::Vector2> g_uvs; // [14] [2]

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


    if (uv_original) {

        uint32_t numUVs = uv_original->size();

        t_group_type = UV_ref->group;

        //////////////////////////////////////////////////////////////////////// adding all uv positions 
        if (uv_original)
        {
            g_uvs.reserve(numUVs);

            for (uint32_t iterUV = 0; iterUV < numUVs; iterUV++) {

                //std::cout << "final UV values at index " << iterUV << ": --> " << uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w << "  " << uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w << std::endl;

                g_uvs.push_back({ uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w, uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w });
            }
        }
    }

    g_nuvs = g_uvs.size();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////////////////////// ADDING ALL VERTEX POSITIONS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    g_verts.reserve(num_vertices);

    for (uint32_t iter = 0; iter < num_vertices; iter++)
    {
        const DD::Image::Vector3& local_point = geoInfo.point_array()[iter];

        DD::Image::Vector4 w_point(geoInfo.matrix * local_point);

        g_verts.push_back({ w_point[0],w_point[1],w_point[2] });

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

    Sdc::SchemeType type = shapescheme;
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(opt.GetVtxBoundaryInterpolation());
    options.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_ALL);
    options.SetCreasingMethod(Sdc::Options::CREASE_CHAIKIN);
    options.SetTriangleSubdivision(Sdc::Options::TRI_SUB_SMOOTH);

    // Instantiate a Far::TopologyRefiner from the descriptor
    Far::TopologyRefiner* refiner =
        Far::TopologyRefinerFactory<Descriptor>::Create(desc,
            Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

    // Uniformly refine the topolgy up to 'maxlevel'
    // note: fullTopologyInLastLevel must be true to work with face-varying data
    {
        Far::TopologyRefiner::UniformOptions refineOptions(maxlevel);
        refineOptions.fullTopologyInLastLevel = true;
        refiner->RefineUniform(refineOptions);
    }

    // Allocate and initialize the 'vertex' primvar data (see tutorial 2 for
    // more details).
    std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal());
    Vertex* verts = &vbuffer[0];
    for (int i = 0; i < g_nverts; ++i) {
        verts[i].SetPosition(g_verts[i][0], g_verts[i][1], g_verts[i][2]);
    }

    // Allocate & initialize the first channel of 'face-varying' primvars (UVs)
    std::vector<FVarVertexUV> fvBufferUV(refiner->GetNumFVarValuesTotal(channelUV));
    FVarVertexUV* fvVertsUV = &fvBufferUV[0];
    for (int i = 0; i < g_nuvs; ++i) {
        fvVertsUV[i].u = g_uvs[i][0];
        fvVertsUV[i].v = g_uvs[i][1];
    }

    //// Allocate & interpolate the 'face-varying' primvar data (colors)
    //std::vector<FVarVertexColor> fvBufferColor(refiner->GetNumFVarValuesTotal(channelColor));
    //FVarVertexColor* fvVertsColor = &fvBufferColor[0];
    //for (int i = 0; i < g_ncolors; ++i) {
    //    fvVertsColor[i].r = g_colors[i][0];
    //    fvVertsColor[i].g = g_colors[i][1];
    //    fvVertsColor[i].b = g_colors[i][2];
    //    fvVertsColor[i].a = g_colors[i][3];
    //}

    // Interpolate both vertex and face-varying primvar data
    Far::PrimvarRefiner primvarRefiner(*refiner);
    Vertex* srcVert = verts;
    FVarVertexUV* srcFVarUV = fvVertsUV;
    //FVarVertexColor* srcFVarColor = fvVertsColor;

    for (int level = 1; level <= maxlevel; ++level) {
        Vertex* dstVert = srcVert + refiner->GetLevel(level - 1).GetNumVertices();
        FVarVertexUV* dstFVarUV = srcFVarUV + refiner->GetLevel(level - 1).GetNumFVarValues(channelUV);
        //FVarVertexColor* dstFVarColor = srcFVarColor + refiner->GetLevel(level - 1).GetNumFVarValues(channelColor);

        primvarRefiner.Interpolate(level, srcVert, dstVert);
        primvarRefiner.LimitFaceVarying(srcFVarUV, dstFVarUV, channelUV);
        primvarRefiner.InterpolateFaceVarying(level, srcFVarUV, dstFVarUV, channelUV);
        //primvarRefiner.InterpolateFaceVarying(level, srcFVarColor, dstFVarColor, channelColor);

        srcVert = dstVert;
        srcFVarUV = dstFVarUV;
        //srcFVarColor = dstFVarColor;
    }

    // Approximate normals
    Far::TopologyLevel const& refLastLevel = refiner->GetLevel(maxlevel);
    int nverts = refLastLevel.GetNumVertices();
    int nfaces = refLastLevel.GetNumFaces();
    int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;

    std::vector<Vertex> normals(nverts);

    // Different ways to approximate smooth normals
    //
    // For details check the description at the beginning of the file
    if (normalApproximation == Limit) {

        // Approximation using the normal at the limit with verts that are 
        // not at the limit
        //
        // For details check the description at the beginning of the file

        std::vector<Vertex> fineLimitPos(nverts);
        std::vector<Vertex> fineDu(nverts);
        std::vector<Vertex> fineDv(nverts);

        primvarRefiner.Limit(&verts[firstOfLastVerts], fineLimitPos, fineDu, fineDv);

        for (int vert = 0; vert < nverts; ++vert) {
            float const* du = fineDu[vert].GetPosition();
            float const* dv = fineDv[vert].GetPosition();

            float norm[3];
            cross(du, dv, norm);
            normals[vert].SetPosition(norm[0], norm[1], norm[2]);
        }

    }
    else if (normalApproximation == CrossQuad) {

        // Approximate smooth normals by accumulating normal vectors computed as
        // the cross product of two vectors generated by the 4 verts that 
        // form each quad
        //
        // For details check the description at the beginning of the file

        for (int f = 0; f < nfaces; f++) {
            Far::ConstIndexArray faceVertices = refLastLevel.GetFaceVertices(f);

            // We will use the first three verts to calculate a normal
            const float* v0 = verts[firstOfLastVerts + faceVertices[0]].GetPosition();
            const float* v1 = verts[firstOfLastVerts + faceVertices[1]].GetPosition();
            const float* v2 = verts[firstOfLastVerts + faceVertices[2]].GetPosition();
            const float* v3 = verts[firstOfLastVerts + faceVertices[3]].GetPosition();

            // Calculate the cross product between the vectors formed by v1-v0 and
            // v2-v0, and then normalize the result
            float normalCalculated[] = { 0.0,0.0,0.0 };
            float a[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };
            float b[3] = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
            cross(a, b, normalCalculated);
            normalize(normalCalculated);

            // Accumulate that normal on all verts that are part of that face
            for (int vInFace = 0; vInFace < faceVertices.size(); vInFace++) {

                int vertexIndex = faceVertices[vInFace];
                normals[vertexIndex].position[0] += normalCalculated[0];
                normals[vertexIndex].position[1] += normalCalculated[1];
                normals[vertexIndex].position[2] += normalCalculated[2];
            }
        }

    }
    else if (normalApproximation == CrossTriangle) {

        // Approximate smooth normals by accumulating normal vectors computed as
        // the cross product of two vectors generated by 3 verts of the quad
        //
        // For details check the description at the beginning of the file

        for (int f = 0; f < nfaces; f++) {
            Far::ConstIndexArray faceVertices = refLastLevel.GetFaceVertices(f);

            // We will use the first three verts to calculate a normal
            const float* v0 = verts[firstOfLastVerts + faceVertices[0]].GetPosition();
            const float* v1 = verts[firstOfLastVerts + faceVertices[1]].GetPosition();
            const float* v2 = verts[firstOfLastVerts + faceVertices[2]].GetPosition();

            // Calculate the cross product between the vectors formed by v1-v0 and
            // v2-v0, and then normalize the result
            float normalCalculated[] = { 0.0,0.0,0.0 };
            float a[3] = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
            float b[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };
            cross(a, b, normalCalculated);
            normalize(normalCalculated);

            // Accumulate that normal on all verts that are part of that face
            for (int vInFace = 0; vInFace < faceVertices.size(); vInFace++) {

                int vertexIndex = faceVertices[vInFace];
                normals[vertexIndex].position[0] += normalCalculated[0];
                normals[vertexIndex].position[1] += normalCalculated[1];
                normals[vertexIndex].position[2] += normalCalculated[2];
            }
        }
    }

    // Finally we just need to normalize the accumulated normals
    for (int vert = 0; vert < nverts; ++vert) {
        normalize(&normals[vert].position[0]);
    }

#ifdef  DEBUG

    { // Output OBJ of the highest level refined -----------

        // Print vertex positions
        for (int vert = 0; vert < nverts; ++vert) {
            float const* pos = verts[firstOfLastVerts + vert].GetPosition();
            printf("v %f %f %f\n", pos[0], pos[1], pos[2]);
        }

        // Print vertex normals
        for (int vert = 0; vert < nverts; ++vert) {
            float const* pos = normals[vert].GetPosition();
            printf("vn %f %f %f\n", pos[0], pos[1], pos[2]);
        }

        // Print uvs
        int nuvs = refLastLevel.GetNumFVarValues(channelUV);
        int firstOfLastUvs = refiner->GetNumFVarValuesTotal(channelUV) - nuvs;
        for (int fvvert = 0; fvvert < nuvs; ++fvvert) {
            FVarVertexUV const& uv = fvVertsUV[firstOfLastUvs + fvvert];
            printf("vt %f %f\n", uv.u, uv.v);
        }

        // Print faces
        for (int face = 0; face < nfaces; ++face) {
            Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(face);
            Far::ConstIndexArray fuvs = refLastLevel.GetFaceFVarValues(face, channelUV);

            // all refined Catmark faces should be quads
            assert(fverts.size() == 4 && fuvs.size() == 4);

            printf("f ");
            for (int vert = 0; vert < fverts.size(); ++vert) {
                // OBJ uses 1-based arrays...
                printf("%d/%d/%d ", fverts[vert] + 1, fuvs[vert] + 1, fverts[vert] + 1);
            }
            printf("\n");
        }
    }

    std::cout << "geometry_list out: \n";
#endif //  DEBUG


    if (op->rebuild(Mask_Primitives)) {
        out.add_object(obj);
        //out.writable_primitive(obj, );

        auto polymesh = new PolyMesh(nfaces, 4);
        for (size_t i = 0; i < nfaces; ++i) {
            Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(i);
            std::vector<int> face;
            face.reserve(4);
            for (int i = 0; i < 4; ++i) face.push_back(fverts[i]);
            polymesh->add_face(4, face.data());
        }
        out.add_primitive(obj, polymesh);
        

        //// adding material
        out[obj].material = geoInfo.material;
        out[obj].useMaterialContext = geoInfo.useMaterialContext;
        out[obj].materialContext = geoInfo.materialContext;

        // Force points and attributes to update:
        op->set_rebuild(Mask_Points | Mask_Attributes);
    }

    if (op->rebuild(Mask_Points)) {
        PointList& points = *out.writable_points(obj);
        points.clear();
        points.resize(nverts);


        for (int i = 0; i < nverts; ++i) {
            float const* pos = verts[firstOfLastVerts + i].GetPosition();
            points[i] = { pos[0], pos[1], pos[2] };
#ifdef  DEBUG
            printf("vertices_out %f %f %f\n", pos[0], pos[1], pos[2]);
#endif
        }
    }

    int nuvs = refLastLevel.GetNumFVarValues(channelUV);
    int firstOfLastUvs = refiner->GetNumFVarValuesTotal(channelUV) - nuvs;
    if (op->rebuild(Mask_Attributes) && nuvs != 0) {
        Attribute* uv = out.writable_attribute(obj, Group_Vertices, "uv", VECTOR4_ATTRIB);
        assert(uv != nullptr);
        uv->clear();
        uv->resize(nuvs * nfaces);

        int faceCount = 0;
        for (int face = 0; face < nfaces; ++face) {
            Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(face);
            Far::ConstIndexArray fuvs = refLastLevel.GetFaceFVarValues(face, channelUV);

            // all refined Catmark faces should be quads
            assert(fverts.size() == 4 && fuvs.size() == 4);



            for (int fvvert = 0; fvvert < fverts.size(); ++fvvert) {
                FVarVertexUV const& uvIn = fvVertsUV[firstOfLastUvs + fuvs[fvvert]];
                uv->vector4(faceCount + fvvert).set(uvIn.u, uvIn.v, 0, 1);
#ifdef  DEBUG
                printf("uv_out %d\n", fuvs[fvvert]);
                printf("uv_out %f %f\n", uvIn.u, uvIn.v);
#endif
            }
            faceCount += 4;
        }
    }

    if (op->rebuild(Mask_Attributes)) {
        Attribute* N = out.writable_attribute(obj, Group_Points, "N", NORMAL_ATTRIB);
        assert(N != nullptr);
        N->clear();
        N->resize(nverts);


        for (int vert = 0; vert < nverts; ++vert) {
            float const* pos = normals[vert].GetPosition();
            N->normal(vert).set(pos[0], pos[1], pos[2]);
        }
    }

    delete refiner;

}


