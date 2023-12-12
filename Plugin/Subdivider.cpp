#include "Subdivider.h"

//  Local headers with support for this tutorial in "namespace tutorial"
//#include "./meshLoader.h"

//#include "./objWriter.h"

Shape* read_GeoOp_input(const DD::Image::GeoInfo& geoInfo, DD::Image::Op* op, Scheme shapescheme, bool isLeftHanded)
{
    Shape* s = new Shape;

    s->scheme = shapescheme;
    s->isLeftHanded = isLeftHanded;

    const DD::Image::PointList* geoPoints = geoInfo.point_list();

    // reserve mesh elements
    uint32_t num_vertices = (uint32_t)geoPoints->size();;

#ifdef DEBUG_MESH
    std::cout << "readgeo --> total vertex amount: " << num_vertices << std::endl;
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------> , UVs, Normals(disabled) vertex
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////get normal information
    const DD::Image::AttribContext* N_ref = geoInfo.get_attribcontext("N");
    const DD::Image::AttributePtr vertex_normal = N_ref ? N_ref->attribute : DD::Image::AttributePtr();

    if (!vertex_normal) {
        op->error("Missing \"N\" channel from geometry");
        return false;
    }

    ////////////////////////////////////////////////////////////////////// adding all normals positions 
    if (vertex_normal) {

        uint32_t numNormals = vertex_normal->size();
        s->normals.reserve(numNormals * 3);

        for (int iterN = 0; iterN < numNormals; iterN++) {

            const DD::Image::Vector3& localN = vertex_normal->normal(iterN);
            DD::Image::Vector4 worldN = geoInfo.matrix * localN;

            //std::cout << "vertex N values at index " << iterN << ": --> " << localN.x << "  " << localN.y << "  " << localN.z << std::endl;

            s->normals.push_back(worldN.x);
            s->normals.push_back(worldN.y);
            s->normals.push_back(worldN.z);

        }
    }
    DD::Image::GroupType n_group_type = DD::Image::Group_None;
    if(vertex_normal)
        DD::Image::GroupType n_group_type = N_ref->group;    // normal group type 

    // get the original uv attribute used to restore untouched uv coordinate
    //const DD::Image::AttribContext* UV_ref = geoInfo.UV_ref;
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
            s->uvs.reserve(numUVs * 2);

            for (uint32_t iterUV = 0; iterUV < numUVs; iterUV++) {

                //std::cout << "final UV values at index " << iterUV << ": --> " << uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w << "  " << uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w << std::endl;

                s->uvs.push_back(uv_original->vector4(iterUV).x / uv_original->vector4(iterUV).w);
                s->uvs.push_back(uv_original->vector4(iterUV).y / uv_original->vector4(iterUV).w);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////////////////////// ADDING ALL VERTEX POSITIONS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    s->verts.reserve(num_vertices * 3);

    for (uint32_t iter = 0; iter < num_vertices; iter++)
    {
        const DD::Image::Vector3& local_point = geoInfo.point_array()[iter];

        DD::Image::Vector4 w_point(geoInfo.matrix * local_point);

        for (uint32_t v = 0; v < 3; ++v)
            s->verts.push_back(w_point[v]);

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////--------------------------------------------------------------------> GENERATE MESH
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///tri or quads primitives
    const DD::Image::Primitive** prim_array = geoInfo.primitive_array();
    if (prim_array)
    {
        const int num_prims = geoInfo.primitives(); //amount of primitives


        //// the following code stores the face numbers that are shared by each vertex
        //vertex_faces.reserve(num_vertices);
        //std::vector<int> faceList;
        //faceList.reserve(num_vertices);

        //for (uint32_t j = 0; j < num_vertices; j++)
        //	vertex_faces.push_back(faceList); // store empty lists so that we can use the [] notation below

        // vertex index array
        std::vector<uint32_t> varray;

        int vertexCount = 0; // faces * vertex ----- for uvs index
        //int count_for_normals = 0; // faces count for normals smoothing

        for (int p = 0; p < num_prims; ++p) {

            // get primitive
            const DD::Image::Primitive* prim = geoInfo.primitive(p);
            const uint32_t num_faces = prim->faces();
            //std::cout << "total sub faces amount: " << num_faces << std::endl;

            // iterate all faces in the primitive
            for (uint32_t f = 0; f < num_faces; f++) {

                // get the number of vertex used by face f
                const uint32_t num_verts = prim->face_vertices(f);

                s->nvertsPerFace.push_back(num_verts);

                // get all used vertices index
                varray.resize(num_verts);

                prim->get_face_vertices(f, &varray[0]);
                //if(n_group_type == DD::Image::Group_Vertices)
                //    std::swap(varray[2], varray[3]);

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
                    unsigned ni = (n_group_type == DD::Image::Group_Points) ? vi : (vertexOffset + v);

                    s->faceverts.push_back(vi);
                    if (uv_original) s->faceuvs.push_back(ti);
                    if (vertex_normal) s->facenormals.push_back(ni);

                }

                vertexCount += num_verts;

            }
        }

    }

    return s;
}

std::vector<ShapeData> tessellateToObj(Far::TopologyRefiner const& meshTopology, std::vector<float> const& meshVtxData, int vtxDataSize, std::vector<float> const& meshFVarData, int fvarDataSize, const bool tessQuadsFlag, const int tessUniformRate) {

    //
    //  Use simpler local type names for the Surface and its factory:
    //
    typedef Bfr::RefinerSurfaceFactory<> SurfaceFactory;
    typedef Bfr::Surface<float>          Surface;
    typedef Surface::PointDescriptor     SurfacePoint;

    //
    //  Identify the source positions and UVs within more general data
    //  arrays for the mesh. If position and/or UV are not at the start
    //  of the vtx and/or fvar data, simply offset the head of the array
    //  here accordingly:
    //
    bool meshHasUVs = (meshTopology.GetNumFVarChannels() > 0);

    float const* meshPosData = meshVtxData.data();
    SurfacePoint  meshPosPoint(3, vtxDataSize);

    float const* meshUVData = meshHasUVs ? meshFVarData.data() : 0;
    SurfacePoint  meshUVPoint(2, fvarDataSize);

    //
    //  Initialize the SurfaceFactory for the given base mesh (very low
    //  cost in terms of both time and space) and tessellate each face
    //  independently (i.e. no shared vertices):
    //
    //  Note that the SurfaceFactory is not thread-safe by default due to
    //  use of an internal cache.  Creating a separate instance of the
    //  SurfaceFactory for each thread is one way to safely parallelize
    //  this loop.  Another (preferred) is to assign a thread-safe cache
    //  to the single instance.
    //
    //  First declare any evaluation options when initializing:
    //
    //  When dealing with face-varying data, an identifier is necessary
    //  when constructing Surfaces in order to distinguish the different
    //  face-varying data channels. To avoid repeatedly specifying that
    //  identifier when only one is present (or of interest), it can be
    //  specified via the Options.
    //
    SurfaceFactory::Options surfaceOptions;
    if (meshHasUVs) {
        surfaceOptions.SetDefaultFVarID(0);
    }

    SurfaceFactory surfaceFactory(meshTopology, surfaceOptions);

    //
    //  The Surface to be constructed and evaluated for each face -- as
    //  well as the intermediate and output data associated with it -- can
    //  be declared in the scope local to each face. But since dynamic
    //  memory is involved with these variables, it is preferred to declare
    //  them outside that loop to preserve and reuse that dynamic memory.
    //
    Surface posSurface;
    Surface uvSurface;

    std::vector<float> facePatchPoints;

    //std::vector<float> outCoords;
    //std::vector<float> outPos, outDu, outDv;
    //std::vector<float> outUV;
    //std::vector<int>   outFacets;

    std::vector<ShapeData> shape_list;

    //
    //  Assign Tessellation Options applied for all faces.  Tessellations
    //  allow the creating of either 3- or 4-sided faces -- both of which
    //  are supported here via a command line option:
    //
    int const tessFacetSize = 3 + tessQuadsFlag;

    Bfr::Tessellation::Options tessOptions;
    tessOptions.SetFacetSize(tessFacetSize);
    tessOptions.PreserveQuads(tessQuadsFlag);

    //
    //  Process each face, writing the output of each in Obj format:
    //
    //tutorial::ObjWriter objWriter(options.outputObjFile);

    int numFaces = surfaceFactory.GetNumFaces();
    for (int faceIndex = 0; faceIndex < numFaces; ++faceIndex) {
        //
        //  Initialize the Surfaces for position and UVs of this face.
        //  There are two ways to do this -- both illustrated here:
        //
        //  Creating Surfaces for the different data interpolation types
        //  independently is clear and convenient, but considerable work
        //  may be duplicated in the construction process in the case of
        //  non-linear face-varying Surfaces. So unless it is known that
        //  face-varying interpolation is linear, use of InitSurfaces()
        //  is generally preferred.
        //
        //  Remember also that the face-varying identifier is omitted from
        //  the initialization methods here as it was previously assigned
        //  to the SurfaceFactory::Options. In the absence of an assignment
        //  of the default FVarID to the Options, a failure to specify the
        //  FVarID here will result in failure.
        //
        //  The cases below are expanded for illustration purposes, and
        //  validity of the resulting Surface is tested here, rather than
        //  the return value of initialization methods.
        //
        bool createSurfacesTogether = true;
        if (!meshHasUVs) {
            surfaceFactory.InitVertexSurface(faceIndex, &posSurface);
        }
        else if (createSurfacesTogether) {
            surfaceFactory.InitSurfaces(faceIndex, &posSurface, &uvSurface);
        }
        else {
            if (surfaceFactory.InitVertexSurface(faceIndex, &posSurface)) {
                surfaceFactory.InitFaceVaryingSurface(faceIndex, &uvSurface);
            }
        }
        if (!posSurface.IsValid()) continue;

        //
        //  Declare a simple uniform Tessellation for the Parameterization
        //  of this face and identify coordinates of the points to evaluate:
        //
        Bfr::Tessellation tessPattern(posSurface.GetParameterization(),
            tessUniformRate, tessOptions);

        int numOutCoords = tessPattern.GetNumCoords();

        ShapeData faceOut;

        faceOut.outCoords.resize(numOutCoords * 2);

        tessPattern.GetCoords(faceOut.outCoords.data());

        //
        //  Prepare the patch points for the Surface, then use them to
        //  evaluate output points for all identified coordinates:
        //
        //  Evaluate vertex positions:
        {
            //  Resize patch point and output arrays:
            int pointSize = meshPosPoint.size;

            facePatchPoints.resize(posSurface.GetNumPatchPoints() * pointSize);

            faceOut.outPos.resize(numOutCoords* pointSize);
            faceOut.outDu.resize(numOutCoords* pointSize);
            faceOut.outDv.resize(numOutCoords* pointSize);

            //  Populate patch point and output arrays:
            float* patchPosData = facePatchPoints.data();
            SurfacePoint  patchPosPoint(pointSize);

            posSurface.PreparePatchPoints(meshPosData, meshPosPoint,
                patchPosData, patchPosPoint);

            for (int i = 0, j = 0; i < numOutCoords; ++i, j += pointSize) {
                posSurface.Evaluate(&faceOut.outCoords[i * 2],
                    patchPosData, patchPosPoint,
                    &faceOut.outPos[j], &faceOut.outDu[j], &faceOut.outDv[j]);
            }
        }

        //  Evaluate face-varying UVs (when present):
        if (meshHasUVs) {
            //  Resize patch point and output arrays:
            //      - note reuse of the same patch point array as position
            int pointSize = meshUVPoint.size;

            facePatchPoints.resize(uvSurface.GetNumPatchPoints() * pointSize);

            faceOut.outUV.resize(numOutCoords * pointSize);

            //  Populate patch point and output arrays:
            float* patchUVData = facePatchPoints.data();
            SurfacePoint patchUVPoint(pointSize);

            uvSurface.PreparePatchPoints(meshUVData, meshUVPoint,
                patchUVData, patchUVPoint);

            for (int i = 0, j = 0; i < numOutCoords; ++i, j += pointSize) {
                uvSurface.Evaluate(&faceOut.outCoords[i * 2],
                    patchUVData, patchUVPoint,
                    &faceOut.outUV[j]);
            }
        }

        //
        //  Identify the faces of the Tessellation:
        //
        //  Note the need to offset vertex indices for the output faces --
        //  using the number of vertices generated prior to this face. One
        //  of several Tessellation methods to transform the facet indices
        //  simply translates all indices by the desired offset.
        //
        int objVertexIndexOffset = 0; // = objWriter.GetNumVertices();
        for (const auto& s : shape_list) 
            objVertexIndexOffset += s.outPos.size() / (3 + tessQuadsFlag);
        

        faceOut.numFacets = tessPattern.GetNumFacets();
        faceOut.outFacets.resize(faceOut.numFacets * tessFacetSize);
        tessPattern.GetFacets(faceOut.outFacets.data());

        tessPattern.TransformFacetCoordIndices(faceOut.outFacets.data(),
            objVertexIndexOffset);

        ////
        ////  Write the evaluated points and faces connecting them as Obj:
        ////
        //objWriter.WriteGroupName("baseFace_", faceIndex);

        //if (meshHasUVs && options.uv2xyzFlag) {
        //    objWriter.WriteVertexPositions(outUV, 2);
        //    objWriter.WriteFaces(outFacets, tessFacetSize, false, false);
        //}
        //else {
        //    objWriter.WriteVertexPositions(outPos);
        //    objWriter.WriteVertexNormals(outDu, outDv);
        //    if (meshHasUVs) {
        //        objWriter.WriteVertexUVs(outUV);
        //    }
        //    objWriter.WriteFaces(outFacets, tessFacetSize, true, meshHasUVs);
        //}

        shape_list.push_back(std::move(faceOut));
    }
    return shape_list;
}

Far::TopologyRefiner* readTopologyRefiner(const DD::Image::GeoInfo& geoInfo, DD::Image::Op* op, Sdc::SchemeType schemeType, std::vector<float>& posVector, std::vector<float>& uvVector)
{

    const Shape* shape = 0;

    //std::ifstream ifs(filename);
    //if (ifs) {
    //    std::stringstream ss;
    //    ss << ifs.rdbuf();
    //    ifs.close();
    //    std::string shapeString = ss.str();

    shape = read_GeoOp_input(geoInfo, op);
    if (shape == NULL) {
        fprintf(stderr, "no geo found");
        return nullptr;
    }

    Sdc::SchemeType sdcType = GetSdcType(*shape);
    Sdc::Options    sdcOptions = GetSdcOptions(*shape);

    Far::TopologyRefiner* refiner = Far::TopologyRefinerFactory<Shape>::Create(
        *shape, Far::TopologyRefinerFactory<Shape>::Options(sdcType, sdcOptions));
    if (refiner == 0) {
        fprintf(stderr,
            "Error:  Unable to construct TopologyRefiner from Obj file '%s'\n",
            "geoInfo");
        return 0;
    }

    int numVertices = refiner->GetNumVerticesTotal();
    posVector.resize(numVertices * 3);
    std::memcpy(&posVector[0], &shape->verts[0], 3 * numVertices * sizeof(float));

    uvVector.resize(0);
    if (refiner->GetNumFVarChannels()) {
        int numUVs = refiner->GetNumFVarValuesTotal(0);
        uvVector.resize(numUVs * 2);
        std::memcpy(&uvVector[0], &shape->uvs[0], 2 * numUVs * sizeof(float));
    }

    delete shape;
    return refiner;
}



