/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team, bachelor student Marcel Heine
 * All rights reserved.
 */

#include "ViewOptimizationRenderer.h"

#include "mmstd_gl/view/View3DGL.h"
#include "protein_calls/MolecularDataCall.h"
#include "compositing_gl/CompositingCalls.h"

using namespace megamol::protein_gl;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getTargetMeshData_("getTargetMeshData",
              "Connects the renderer to a data provider to retrieve target mesh data")
        , getLigandPDBData_("getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , getTexture_("InputTexture", "Access texture that is used to calulate the Viewpoint Entropy")
        , _cutTriangleMesh("cutTriangleMesh", "Forwards either the mesh data of a previous 'MSMSMeshLoader' or a reduced version")
        , optimizeCamera("optimizeCamera", "Acts as a button to set the camera to view the ligand binding site")
        , refreshTargetMesh("refreshTargetMesh", "Acts as a button to refresh/recalculate the rendered target triangle mesh")
        , renderTargetMeshModeParam("targetMeshRenderMode", "The target mesh rendering mode.")
        , molRadiusSummand("molRadiusSummand", "Value to be added to the ligand molecule radius. Represents atom radius and target-ligand gap distance.")
        , currentTargetMeshData()
        , currentTargetMeshRenderMode(WHOLE_MESH)
        , bbox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f)
        , datahash(0)
        , reload_colors(true) {
    // parameters
    this->optimizeCamera.SetParameter(new core::param::BoolParam(false));
    this->MakeSlotAvailable(&this->optimizeCamera);

    this->refreshTargetMesh.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->refreshTargetMesh);

    megamol::core::param::EnumParam* ctmrm = new megamol::core::param::EnumParam(int(this->currentTargetMeshRenderMode));
    ctmrm->SetTypePair(WHOLE_MESH, "whole mesh");
    ctmrm->SetTypePair(NAIVE_CAVETY, "naive cavety");
    this->renderTargetMeshModeParam << ctmrm;
    this->MakeSlotAvailable(&this->renderTargetMeshModeParam);

    this->molRadiusSummand.SetParameter(new core::param::FloatParam(2.0f));  // value approximates the average gap between ligand and target + atom radius. (arbitrarily chosen)
    this->MakeSlotAvailable(&this->molRadiusSummand);

    // data in slot(s)
    this->getTargetMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    this->getTargetMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&getTargetMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
    this->getLigandPDBData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->getLigandPDBData_);

    this->getTexture_.SetCompatibleCall<megamol::compositing_gl::CallTexture2DDescription>();
    this->MakeSlotAvailable(&this->getTexture_);

    // data out slot(s)
    this->_cutTriangleMesh.SetCallback(geocalls_gl::CallTriMeshDataGL::ClassName(), "GetData", &ViewOptimizationRenderer::getDataCallback);
    this->_cutTriangleMesh.SetCallback(geocalls_gl::CallTriMeshDataGL::ClassName(), "GetExtent", &ViewOptimizationRenderer::getExtentCallback);
    this->MakeSlotAvailable(&this->_cutTriangleMesh);
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
    if (this->optimizeCamera.Param<core::param::BoolParam>()->Value()) {

        // Calculate naive camera position based on ligand center 
        /* ------- Get Ligand Data ------- */
        protein_calls::MolecularDataCall* ligand = this->getLigandPDBData_.CallAs<protein_calls::MolecularDataCall>();
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

        float radius = moleculeRadius(ligand, pos0, ligandCenter, this->molRadiusSummand.Param<core::param::FloatParam>()->Value());

        /* ------- Get Naive Camera Direction ------- */

        // get the naive camera direction
        glm::vec3 naiveCamDirection = naiveCameraDirection(target, ligandCenter, radius);

        /* ------- Viewpoint Entropy: Import Rendered Texture ------- */

        std::shared_ptr<glowl::Texture2D> colorTexture = nullptr;
        megamol::compositing_gl::CallTexture2D* ct2d = this->getTexture_.CallAs<megamol::compositing_gl::CallTexture2D>();
        if (ct2d != NULL) {
            (*ct2d)(0);
            colorTexture = ct2d->getData();
        }
        if (colorTexture == nullptr) {
            return false;
        }

        const unsigned int hight = colorTexture->getHeight();
        const unsigned int width = colorTexture->getWidth();
        const unsigned int rgb = 3;
        const unsigned int dataTypeBits = 1;
        const unsigned int storageSize = hight * width * rgb * dataTypeBits;
        uint8_t* textureData = new uint8_t[storageSize];

        glGetTextureImage(colorTexture->getName(), 0, GL_RGB, GL_UNSIGNED_BYTE, storageSize, textureData);

        for (unsigned int i = storageSize/6; i < storageSize/6 + 16; i++) {
            //std::cout << "Texture data" << i << " : " << int(textureData[i]) << "\n";
            std::cout << "Texture data center" << i + storageSize/3 << " : " << unsigned int(textureData[i+storageSize / 3]) << "\n";
            //std::cout << "Texture data" << i + storageSize / 3 * 2 << " : " << int(textureData[i + storageSize / 3 * 2]) << "\n";
        }

        for (unsigned int i = 0; i < 9; i++) {
            std::cout << "Texture data begining" << i << " : " << unsigned int(textureData[i]) << "\n";
        }

        for (unsigned int i = storageSize - 15; i < storageSize; i++) {
            std::cout << "Texture data end" << i << " : " << unsigned int(textureData[i]) << "\n";
        }

        std::cout << "Backgroundcolor: "
            << call.BackgroundColor().r << ", "
            << call.BackgroundColor().g << ", "
            << call.BackgroundColor().b << ", "
            << call.BackgroundColor().a;

        /* ------- Set New Camera ------- */

        const float camDistFactor = 16.0f;
        
        newCallCamera(call, naiveCamDirection, ligandCenter, radius, camDistFactor);

        /* ------- Cleanup & Global Value(s) ------- */

        // delete pointers and ensure, that this if branch is executed once until user request
        this->optimizeCamera.Param<core::param::BoolParam>()->SetValue(false);     
        delete[] pos0, textureData;
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

/* =============== Necessary Methods For Target Data Relay =============== */

bool ViewOptimizationRenderer::getDataCallback(core::Call& caller) {
    /* ------- Get Call Data / Prepare Calls ------- */
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

    // Update 'currentTargetMeshRenderMode' if necessary 
    this->UpdateParameters();

    /* ------- Mesh Rendering ------- */

    // initialize or reload the container 'currentTargetMeshData' with target mesh data
    if (this->refreshTargetMesh.Param<core::param::BoolParam>()->Value()) {
        // decide which mesh will be passed on to rendering
        if (this->currentTargetMeshRenderMode == NAIVE_CAVETY) {
            protein_calls::MolecularDataCall* ligand = this->getLigandPDBData_.CallAs<protein_calls::MolecularDataCall>();
            if (ligand == NULL)
                return false;

            if (!(*ligand)(protein_calls::MolecularDataCall::CallForGetData))
                return false;
            if (ligand->AtomCount() == 0)
                return true;

            const unsigned int ligandAtomCount = ligand->AtomCount();
            float* pos0 = new float[ligandAtomCount * 3];
            memcpy(pos0, ligand->AtomPositions(), ligandAtomCount * 3 * sizeof(float));

            glm::vec3 ligandCenter = moleculeCenter(pos0, ligandAtomCount);
            float radius = moleculeRadius(ligand, pos0, ligandCenter, this->molRadiusSummand.Param<core::param::FloatParam>()->Value());

            this->currentTargetMeshData = naiveCavetyCutter(targetCtmd->Objects()[0], ligandCenter, radius, true);
            delete[] pos0;
        }
        else if (this->currentTargetMeshRenderMode == WHOLE_MESH) {
            this->currentTargetMeshData = targetCtmd->Objects();
        }     

        this->refreshTargetMesh.Param<core::param::BoolParam>()->SetValue(false);
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

/* =============== Parameters =============== */

/*
 * Updates all module parameters 
 */
void ViewOptimizationRenderer::UpdateParameters() {
    // target mesh rendering mode param
    static MeshRenderMode lastMode = WHOLE_MESH;        // Might cause problems when multiple instances of 'ViewOptimizationRenderer' are run simultaneously
    if (this->renderTargetMeshModeParam.IsDirty()) {
        this->currentTargetMeshRenderMode =
            static_cast<MeshRenderMode>(int(this->renderTargetMeshModeParam.Param<megamol::core::param::EnumParam>()->Value()));

        // Ensure, that the mesh is only updated when requested 
        if (this->currentTargetMeshRenderMode != lastMode) {
            this->refreshTargetMesh.Param<core::param::BoolParam>()->SetValue(true);
            lastMode = this->currentTargetMeshRenderMode;
        }
    }
}

/* =============== Supporting Functions =============== */

glm::vec3 ViewOptimizationRenderer::moleculeCenter(float* positions, unsigned int atomCount) {
    glm::vec3 ligandCenter = glm::vec3(0, 0, 0);
    for (unsigned int i = 0; i < atomCount * 3; i++) {
        ligandCenter[i % 3] += positions[i] / atomCount;
    }

    return ligandCenter;
}

float ViewOptimizationRenderer::moleculeRadius(megamol::protein_calls::MolecularDataCall* mol, float* atomPos, glm::vec3 ligandCenter, float buffer) {
    // get radius, which will be used to select vertices around the ligand center
    // i.e get max atom distance from ligand center
    float radius = 0;
    for (int i = 0; i < mol->AtomCount(); i++) {
        float dist = glm::distance(glm::vec3(atomPos[i * 3], atomPos[i * 3 + 1], atomPos[i * 3 + 2]), ligandCenter);
        if (dist > radius) {
            radius = dist;
        }
    }
    radius += buffer;

    return radius;
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
    core::view::Camera::Pose newCamPose = cam.getPose();  // old camera as default

    // camera-ligandcenter distance is 'camDistFactor'-times the vertex selection radius
    newCamPose.position = ligCenter + direction * glm::vec3(radius * camDistFactor);
    newCamPose.direction = -direction;

    const auto& camIntrinsics = cam.get<core::view::Camera::PerspectiveParameters>();

    call.SetCamera(core::view::Camera(newCamPose, camIntrinsics));
}

megamol::geocalls_gl::CallTriMeshDataGL::Mesh* ViewOptimizationRenderer::naiveCavetyCutter(
    geocalls_gl::CallTriMeshDataGL::Mesh mesh, glm::vec3 ligCenter, float radius, bool altColAndMesh) {

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
    // assuming 'mesh.GetVertexAttribCount' == 2 with 'AttribID 0' = atomIndex and 'AttribID 1' = values
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
                colours[(i - offset) * 3 + j] = mesh.GetColourPointerByte()[i * 3 + j];   // orriginal colors
                //colours[(i - offset) * 3 + j] = char((j == 0) ? (((i - offset) % 235) + 20) : (j == 1) ? (int((i - offset) / 235) % 235 + 20) : (int((i - offset) / (235 * 235)) % 235 + 20));
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
    const unsigned int* facesOld = mesh.GetTriIndexPointerUInt32();
    std::list<unsigned int> facesList;

    unsigned int vertIndx[3];
    for (unsigned int i = 0; i < faceCount; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            vertIndx[j] = inArray(oldVertIndices, facesOld[i * 3 + j], newVertCount);
        }

        if (vertIndx[0] != newVertCount
            && vertIndx[1] != newVertCount
            && vertIndx[2] != newVertCount) {
            for (unsigned int j = 0; j < 3; j++) {
                facesList.push_back(vertIndx[j]);
            }
            newFaceCount++;
        }
    }
    unsigned int* faces = new unsigned int[newFaceCount * 3];
    std::copy(facesList.begin(), facesList.end(), faces);

    // Redo the arrays, so that every face now has it own 3 vertices 
    if (altColAndMesh) {
        float* verticesAlt = new float[newFaceCount * 9];
        float* normalsAlt = new float[newFaceCount * 9];
        unsigned char* coloursAlt = new unsigned char[newFaceCount * 9];

        unsigned int* atomIndexAlt = new unsigned int[newFaceCount * 3];
        float* valuesAlt = new float[newFaceCount * 3];
        unsigned int* facesAlt = new unsigned int[newFaceCount * 3];

        for (unsigned int i = 0; i < newFaceCount; i++) {   // i faces
            for (unsigned int j = 0; j < 3; j++) {          // j vertices
                for (unsigned int k = 0; k < 3; k++) {      // k koordinates
                    verticesAlt[i * 9 + j * 3 + k] = vertices[faces[i * 3 + j] * 3 + k];
                    normalsAlt[i * 9 + j * 3 + k] = normals[faces[i * 3 + j] * 3 + k];
                    //coloursAlt[i * 9 + j * 3 + k] = (i % 2 == 0) ? ((k == 0) ? 210 : ((k == 1) ? 210 : 50)) : 10;
                    //coloursAlt[i * 9 + j * 3 + k] = colours[faces[i * 3 + j] * 3 + k];
                    coloursAlt[i * 9 + j * 3 + k] = char((k == 0) ? ((i % 235) + 20) :
                        (k == 1) ? (int(i / 235) % 235 + 20) :
                        (int(i / (235 * 235)) % 235 + 20));
                }

                atomIndexAlt[i * 3 + j] = atomIndex[faces[i * 3 + j]];
                valuesAlt[i * 3 + j] = values[faces[i * 3 + j]];
                facesAlt[i * 3 + j] = i * 3 + j;
            }
        }

        // initialize the new mesh with the cut data
        cutMesh->SetVertexData(newFaceCount * 3, verticesAlt, normalsAlt, coloursAlt, NULL, true);
        cutMesh->SetTriangleData(newFaceCount, facesAlt, true);
        cutMesh->AddVertexAttribPointer(atomIndexAlt);
        cutMesh->AddVertexAttribPointer(valuesAlt);
        cutMesh->SetMaterial(NULL);

        delete[] vertices, normals, colours, atomIndex, values, faces;
    }
    else {
        // initialize the new mesh with the cut data
        cutMesh->SetVertexData(newVertCount, vertices, normals, colours, NULL, true);
        cutMesh->SetTriangleData(newFaceCount, faces, true);
        cutMesh->AddVertexAttribPointer(atomIndex);
        cutMesh->AddVertexAttribPointer(values);
        cutMesh->SetMaterial(NULL);
    }

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
