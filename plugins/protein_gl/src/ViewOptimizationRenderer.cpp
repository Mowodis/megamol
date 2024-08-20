/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team, bachelor student Marcel Heine
 * All rights reserved.
 */

#include "ViewOptimizationRenderer.h"

#include <glm/glm.hpp>

#include "geometry_calls_gl/CallTriMeshDataGL.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/view/Camera.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"

#include "mmstd_gl/view/View3DGL.h"
#include "protein_calls/MolecularDataCall.h"

using namespace megamol;
using namespace megamol::protein_gl;
using namespace megamol::protein_calls;

using namespace megamol::core::view;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getTargetMeshData_("getTargetMeshData",
              "Connects the renderer to a data provider to retrieve target mesh data")
        , getLigandPDBData_("getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , _cutTriangleMesh("cutTriangleMesh", "Forwards either the mesh data of a previous 'MSMSMeshLoader' or a reduced version")
        , optimizeCamera_("OptimizeCamera", "Acts as a button to set the camera to view the ligand binding site")
        , currentTargetMeshData()
        , bbox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f)
        , datahash(0)
        , reload_colors(true) {

    this->getTargetMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    this->getTargetMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&getTargetMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
    this->getLigandPDBData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->getLigandPDBData_);

    // the data out slot
    this->_cutTriangleMesh.SetCallback(geocalls_gl::CallTriMeshDataGL::ClassName(), "GetData", &ViewOptimizationRenderer::getDataCallback);
    this->_cutTriangleMesh.SetCallback(geocalls_gl::CallTriMeshDataGL::ClassName(), "GetExtent", &ViewOptimizationRenderer::getExtentCallback);
    this->MakeSlotAvailable(&this->_cutTriangleMesh);

    this->optimizeCamera_.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->optimizeCamera_);
}

ViewOptimizationRenderer::~ViewOptimizationRenderer() {
    this->Release();
}

bool ViewOptimizationRenderer::create() {
    return true;
}

void ViewOptimizationRenderer::release() {}

bool ViewOptimizationRenderer::Render(mmstd_gl::CallRender3DGL& call) {

    /* ============= Cut provide Target Triangel Mesh Data  ============= */

    //if (this->coloringModeParam0.IsDirty() || this->coloringModeParam1.IsDirty() ||
    //    this->colorWeightParam.IsDirty() || this->minGradColorParam.IsDirty() ||
    //    this->midGradColorParam.IsDirty() || this->maxGradColorParam.IsDirty() ||
    //    this->prevTime != int(ctmd->FrameID()) || reload_colors) {

    /* ============= Set camera position ============= */

    // execute this branch for only one tick, then waite until the user sets the parameter to true 
    if (this->optimizeCamera_.Param<core::param::BoolParam>()->Value()) {

        // Calculate naive camera position based on ligand center 
        /* ------- Get Ligand Data ------- */
        MolecularDataCall* ligand = this->getLigandPDBData_.CallAs<MolecularDataCall>();
        if (ligand == NULL)
            return false;

        if (!(*ligand)(protein_calls::MolecularDataCall::CallForGetData))
            return false;
        if (ligand->AtomCount() == 0)
            return true;

        /* ------- Get Target Data ------- */
        geocalls_gl::CallTriMeshDataGL* target =
            this->getTargetMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();
        if (target == nullptr) {
            return false;
        }
        target->SetFrameID(static_cast<int>(call.Time()));
        if (!(*target)(1)) {
            return false;
        }

        target->SetFrameID(static_cast<int>(call.Time()));
        if (!(*target)(0)) {
            return false;
        }

        /* ------- Approximate Ligand Radius ------- */
      
        const unsigned int ligandAtomCount = ligand->AtomCount();
        float* pos0 = new float[ligandAtomCount * 3];
        memcpy(pos0, ligand->AtomPositions(), ligandAtomCount * 3 * sizeof(float));

        glm::vec3 ligandCenter = moleculeCenter(pos0, ligandAtomCount);

        // get radius, which will be used to select vertices around the ligand center
        // i.e get max atom distance from ligand center + atom sphere radius + gap of target-ligand
        float radius = 0;
        for (int i = 0; i < ligand->AtomCount(); i++) {
            float dist = glm::distance(glm::vec3(pos0[i * 3], pos0[i * 3 + 1], pos0[i * 3 + 2]), ligandCenter);
            if (dist > radius) {
                radius = dist;
            }
        }
        // adding value approximatly of average gap between ligand and target (arbitrarily chosen)
        // changes in this value might/should/will significantly change the camer orientation
        const float atomRadiusAndGap = 2.0f;
        radius += atomRadiusAndGap;


        /* ------- Get Naive Camera Direction ------- */

        // get the naive camera direction
        glm::vec3 naiveCamDirection = naiveCameraDirection(target, ligandCenter, radius);

        /* ------- Cut The Target Mesh ------- */
        this->currentTargetMeshData = naiveCavetyCutter(target->Objects()[0], ligandCenter, radius);
        //this->currentTargetMeshData = target->Objects();

        std::cout << "Object atrebute nr: " << currentTargetMeshData->GetVertexAttribCount() << "\n";
        std::cout << "Object triangel Pointers: " << currentTargetMeshData->GetTriIndexPointerUInt32()[4] << "\n";


        std::cout << "Vertex Attribute Count: " << currentTargetMeshData->GetVertexAttribCount() << "\n";
        std::cout << "Vertex Attribute dt 0 : " << currentTargetMeshData->GetVertexAttribDataType(0) << "\n";     // 3 -> DT_UINT32
        std::cout << "Vertex Attribute dt 1 : " << currentTargetMeshData->GetVertexAttribDataType(1) << "\n";     // 6 -> DT_float
        std::cout << "Vertex Attribute 0 : " << currentTargetMeshData->GetVertexAttribPointerUInt32(0)[0] << "\n";
        std::cout << "Vertex Attribute 1 : " << currentTargetMeshData->GetVertexAttribPointerFloat(1)[0] << "\n";

        /* ------- Set New Camera ------- */

        const float camDistFactor = 16.0f;
        
        newCallCamera(call, naiveCamDirection, ligandCenter, radius, camDistFactor);

        /* ------- Cleanup & Global Value(s) ------- */

        // delete pointers and ensure, that this if branch is executed once until user request
        this->optimizeCamera_.Param<core::param::BoolParam>()->SetValue(false);     
        delete[] pos0;
    }

    return true;
}

bool ViewOptimizationRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    auto ctmd = this->getTargetMeshData_.CallAs<geocalls_gl::CallTriMeshDataGL>();
    if (ctmd == nullptr) {
        return false;
    }
    ctmd->SetFrameID(static_cast<int>(call.Time()));
    if (!(*ctmd)(1)) {
        return false;
    }

    call.SetTimeFramesCount(ctmd->FrameCount());
    call.AccessBoundingBoxes().Clear();
    call.AccessBoundingBoxes() = ctmd->AccessBoundingBoxes();

    return true;
}

/* =============== Necessary Methods For Data Relay=============== */

bool ViewOptimizationRenderer::getDataCallback(core::Call& caller) {
    geocalls_gl::CallTriMeshDataGL* ctmd = dynamic_cast<geocalls_gl::CallTriMeshDataGL*>(&caller);
    if (ctmd == NULL)
        return false;

    geocalls_gl::CallTriMeshDataGL* targetCtmd =
        this->getTargetMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();
    if (targetCtmd == nullptr) {
        return false;
    }
    targetCtmd->SetFrameID(float(ctmd->FrameID()));
    if (!(*targetCtmd)(1)) {
        return false;
    }
    targetCtmd->SetFrameID(float(ctmd->FrameID()));
    if (!(*targetCtmd)(0)) {
        return false;
    }

    // initialize the container 'currentTargetMeshData' with the unaltered target mesh data
    static bool isInitialized = false;
    if (!isInitialized) {
        this->currentTargetMeshData = targetCtmd->Objects();
        isInitialized = true;
    }

    // pass on the (possibly cut) target mesh data  
    ctmd->SetObjects(1, this->currentTargetMeshData);

    return true;
}

bool ViewOptimizationRenderer::getExtentCallback(core::Call& caller) {
    geocalls_gl::CallTriMeshDataGL* ctmd = dynamic_cast<geocalls_gl::CallTriMeshDataGL*>(&caller);
    if (ctmd == NULL)
        return false;

    unsigned int frameCnt = 1;

    ctmd->SetDataHash(this->datahash);
    bool refresh = false;
    ctmd->SetExtent(frameCnt, this->bbox.Left(), this->bbox.Bottom(), this->bbox.Back(), this->bbox.Right(),
        this->bbox.Top(), this->bbox.Front());
    
    return true;
}

/* =============== Supporting Functions =============== */

glm::vec3 ViewOptimizationRenderer::moleculeCenter(float* positions, unsigned int atomCount) {
    glm::vec3 ligandCenter = glm::vec3(0, 0, 0);
    for (unsigned int i = 0; i < atomCount * 3; i++) {
        ligandCenter[i % 3] += positions[i] / atomCount;
    }

    return ligandCenter;
}

glm::vec3 ViewOptimizationRenderer::naiveCameraDirection(geocalls_gl::CallTriMeshDataGL* ctmd, const glm::vec3 center, const float radius) {
    // get obj: the (first) object of target caller slot. Used for the vertex data
    const auto& target = ctmd->Objects()[0];
    auto vertCount = target.GetVertexCount();
    auto triCount = target.GetTriCount();

    // find all target tie mesh vertices that are within the 'radius' aroud the ligand center
    glm::vec3 naiveCamDirection = glm::vec3(0, 0, 0);
    unsigned int nrSelectedVertices = 0;

    for (int i = 0; i < vertCount; i++) {
        glm::vec3 currentVertex =
            glm::vec3(target.GetVertexPointerFloat()[i * 3],
                target.GetVertexPointerFloat()[i * 3 + 1],
                target.GetVertexPointerFloat()[i * 3 + 2]);

        if (glm::distance(currentVertex, center) <= radius) {
            naiveCamDirection += glm::normalize(center - currentVertex);
            nrSelectedVertices++;
        }
    }

    naiveCamDirection = glm::normalize(naiveCamDirection);

    return naiveCamDirection;
}

void ViewOptimizationRenderer::newCallCamera(
    mmstd_gl::CallRender3DGL& call, glm::vec3 direction, glm::vec3 ligCenter, float radius, float camDistFactor) {
    // Change the camera position coordinates of the CallRenderer3DGL call
    auto& cam = call.GetCamera();
    Camera::Pose newCamPose = cam.getPose();  // old camera as default

    // camera-ligandcenter distance is 'camDistFactor'-times the vertex selection radius
    newCamPose.position = ligCenter + direction * glm::vec3(radius * camDistFactor);
    newCamPose.direction = -direction;

    const auto& camIntrinsics = cam.get<Camera::PerspectiveParameters>();

    call.SetCamera(Camera(newCamPose, camIntrinsics));
}

geocalls_gl::CallTriMeshDataGL::Mesh* ViewOptimizationRenderer::naiveCavetyCutter(
    geocalls_gl::CallTriMeshDataGL::Mesh mesh, glm::vec3 ligCenter, float radius) {

    geocalls_gl::CallTriMeshDataGL::Mesh* cutMesh = new geocalls_gl::CallTriMeshDataGL::Mesh();

    const unsigned int vertCount = mesh.GetVertexCount();
    unsigned int newVertCount = 0;

    // first: count the number of vertices that lay within the radius around the ligand center
    glm::vec3 currentVertex;
    for (unsigned int i = 0; i < vertCount; i++) {
        currentVertex =
            glm::vec3(mesh.GetVertexPointerFloat()[i * 3],
                mesh.GetVertexPointerFloat()[i * 3 + 1],
                mesh.GetVertexPointerFloat()[i * 3 + 2]);

        if (glm::distance(currentVertex, ligCenter) <= radius) {
            newVertCount++;
        }
    }

    float* vertices = new float[newVertCount * 3];
    float* normals = new float[newVertCount * 3];
    unsigned char* colours = new unsigned char[newVertCount * 3];
    // assumes 'mesh.GetVertexAttribCount' == 2 with 'AttribID 0' = atomIndex and 'AttribID 1 = values'
    unsigned int* atomIndex = new unsigned int[newVertCount];
    float* values = new float[newVertCount];

    unsigned int* oldVertIndices = new unsigned int[newVertCount];

    // second: copy all the used vertices with normals and color to new arrays.
    // note the old vertex indices
    unsigned int offset = 0;
    for (unsigned int i = 0; i < vertCount; i++) {
        currentVertex =
            glm::vec3(mesh.GetVertexPointerFloat()[i * 3],
                mesh.GetVertexPointerFloat()[i * 3 + 1],
                mesh.GetVertexPointerFloat()[i * 3 + 2]);

        if (glm::distance(currentVertex, ligCenter) <= radius) {
            for (unsigned int j = 0; j < 3; j++) {
                vertices[(i - offset) * 3 + j] = mesh.GetVertexPointerFloat()[i * 3 + j];
                normals[(i - offset) * 3 + j] = mesh.GetNormalPointerFloat()[i * 3 + j];
                colours[(i - offset) * 3 + j] = mesh.GetColourPointerByte()[i * 3 + j];
            }

            newVertCount++;
            oldVertIndices[i - offset] = i;
            atomIndex[i - offset] = mesh.GetVertexAttribPointerUInt32(0)[i];
            values[i - offset] = mesh.GetVertexAttribPointerFloat(1)[i];
        }
        else {
            offset++;
        }
    }

    const unsigned int faceCount = mesh.GetTriCount();
    unsigned int newFaceCount = 0;

    // determine which triangles are still valid under the new set of vertices and copy those that are
    const unsigned int* faces = mesh.GetTriIndexPointerUInt32();
    std::list<unsigned int> newFacesList;

    unsigned int vertIndx[3];
    for (unsigned int i = 0; i < faceCount; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            vertIndx[j] = inArray(oldVertIndices, faces[i * 3 + j], newVertCount);
        }

        if (vertIndx[0] != newVertCount
            && vertIndx[1] != newVertCount
            && vertIndx[2] != newVertCount) {
            for (unsigned int j = 0; j < 3; j++) {
                newFacesList.push_back(vertIndx[j]);
            }
            newFaceCount++;
        }
    }
    unsigned int* newFaces = new unsigned int[newFaceCount * 3];
    std::copy(newFacesList.begin(), newFacesList.end(), newFaces);

    // initialize the new mesh with the cut data
    cutMesh->SetVertexData(newVertCount, vertices, normals, colours, NULL, true);
    cutMesh->SetTriangleData(newFaceCount, newFaces, true);
    cutMesh->AddVertexAttribPointer(atomIndex);
    cutMesh->AddVertexAttribPointer(values);
    cutMesh->SetMaterial(NULL);

    delete[] oldVertIndices, faces;

    return cutMesh;
}

unsigned int ViewOptimizationRenderer::inArray(unsigned int* arr, unsigned int element, unsigned int arrSize) {
    for (unsigned int i = 0; i < arrSize; i++) {
        if (arr[i] == element) {
            return i;
        }
    }

    return arrSize;
}
