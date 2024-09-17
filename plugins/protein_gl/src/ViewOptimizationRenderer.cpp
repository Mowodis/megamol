/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team
 * Author: Marcel Heine (Bachelor Student)
 * All rights reserved.
 */

#include "ViewOptimizationRenderer.h"

#include "mmstd_gl/view/View3DGL.h"
#include "protein_calls/MolecularDataCall.h"
#include "compositing_gl/CompositingCalls.h"
#include "protein/Icosphere.h"

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
        , molRadiusSummand("molRadiusSummand", "Distance to be added to the ligand molecule radius given in Angström. Represents atom radius and target-ligand gap distance.")
        , camDistFactor("camDistFactor", "Multiplied with the ligand Radius (Distance center to furthest Atom) to form the camera distance to the ligand center.")
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

    this->camDistFactor.SetParameter(new core::param::FloatParam(12.0f));  // comfortable distance (arbitrarily chosen)
    this->MakeSlotAvailable(&this->camDistFactor);

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

    /* ============= Set Dynamic Data Reloading ============= */

    //if (this->coloringModeParam0.IsDirty() || this->coloringModeParam1.IsDirty() ||
    //    this->colorWeightParam.IsDirty() || this->minGradColorParam.IsDirty() ||
    //    this->midGradColorParam.IsDirty() || this->maxGradColorParam.IsDirty() ||
    //    this->prevTime != int(ctmd->FrameID()) || reload_colors) {
    //if (this->numFramesSlot.IsDirty() || this->numSpheresSlot.IsDirty()) {
    //    this->resetFrameCache();
    //    AnimDataModule::setFrameCount(frameCount);
    //    AnimDataModule::initFrameCache(frameCount);
    //    this->numFramesSlot.ResetDirty();
    //    this->numSpheresSlot.ResetDirty();
    //}

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

        const unsigned int height = colorTexture->getHeight();
        const unsigned int width = colorTexture->getWidth();
        const unsigned int rgb = 3;
        const unsigned int storageSize = height * width * rgb;
        uint8_t* textureData = new uint8_t[storageSize];

        glGetTextureImage(colorTexture->getName(), 0, GL_RGB, GL_UNSIGNED_BYTE, storageSize, textureData);

        /* ------- Viewpoint Entropy: Calculation ------- */

        unsigned int* pixelForEachTriangle = pixelColorCounter(textureData, currentTargetMeshData->GetTriCount(), storageSize);      // output will be + 1, due to the background

        float a_t = height * width;      // since the background also counts as a triangle, otherwise sum up visible triangel surface
        float a_i = 0;
        float viewpointEntropy = 0;

        if (a_t != 0) {
            for (unsigned int triangle = 0; triangle < currentTargetMeshData->GetTriCount() + 1; triangle++) {
                a_i = pixelForEachTriangle[triangle];

                if (a_i != 0) {
                    viewpointEntropy += -a_i / a_t * log2(a_i / a_t);
                }
            }
        }


        Icosphere sampleSphere = Icosphere(1.0f, 2, false);

        //const float* ssVertices = sampleSphere.getVertices();
        const float* ssVertices = sampleSphere.getVertices();

        std::cout << "Sample Sphere vertices: " << "\n";
        for (int i = 0; i < sampleSphere.getVertexCount(); i++) {
            std::cout << "Vertex " << i << " : " << ssVertices[i * 3] << ", " << ssVertices[i * 3 + 1] << ", " << ssVertices[i * 3 + 2] << "\n";
        }

        std::cout << "Subdivision: " << sampleSphere.getSubdivision() << "\n";
        unsigned int truevertexCount = pow(sampleSphere.getSubdivision(), 2) * 10 + 2;

        std::cout << "Vertex Count: " << truevertexCount << "\n";
        float* sssVertices = removeDuplicatVertices(ssVertices, sampleSphere.getVertexCount(), sampleSphere.getVertexCount());

        std::cout << "Sample Sphere vertices: " << "\n";
        for (int i = 0; i < sampleSphere.getVertexCount(); i++) {
            std::cout << "Vertex " << i << " : " << sssVertices[i*3] << ", " << sssVertices[i*3+1] << ", " << sssVertices[i*3+2] << "\n";
        }

        std::cout << "a_t " << a_t << "\n";
        std::cout << "Viewpoint Entropy: " << viewpointEntropy << "\n";

        //const glm::vec3* tempDirection = new glm::vec3[1]{ naiveCamDirection };
        //float* tempViewEntropResultes = evaluateViewpoints(call, tempDirection, 1, ligandCenter, radius);

        //std::cout << "VE test: " << tempViewEntropResultes[0] << "\n";
        //delete[] tempDirection, tempViewEntropResultes;


        /* ------- Set New Camera ------- */
        
        setCallCamera(call, naiveCamDirection, ligandCenter, radius, this->camDistFactor.Param<core::param::FloatParam>()->Value());

        /* ------- Cleanup & Global Value(s) ------- */

        // delete pointers and ensure, that this if branch is executed once until user request
        this->optimizeCamera.Param<core::param::BoolParam>()->SetValue(false);     
        delete[] pos0, textureData, ssVertices, sssVertices;
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
    if (this->refreshTargetMesh.Param<core::param::BoolParam>()->Value() || currentTargetMeshData == nullptr) {
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

void ViewOptimizationRenderer::setCallCamera(
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

    float* vertices = new float[newVertCount * 3]{ 0 };
    float* normals = new float[newVertCount * 3]{ 0 };
    unsigned char* colours = new unsigned char[newVertCount * 3]{ 0 };
    // assuming 'mesh.GetVertexAttribCount' == 2 with 'AttribID 0' = atomIndex and 'AttribID 1' = values
    unsigned int* atomIndex = new unsigned int[newVertCount]{ 0 };
    float* values = new float[newVertCount]{ 0 };

    unsigned int* oldVertIndices = new unsigned int[newVertCount]{ 0 };

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

    const unsigned int triangleCount = mesh.GetTriCount();
    unsigned int newTriangleCount = 0;

    // determine which triangles are still valid under the new set of vertices and copy those that are
    const unsigned int* trianglesOld = mesh.GetTriIndexPointerUInt32();
    std::list<unsigned int> trianglesList;

    unsigned int vertIndx[3];
    for (unsigned int i = 0; i < triangleCount; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            vertIndx[j] = inArray(oldVertIndices, trianglesOld[i * 3 + j], newVertCount);
        }

        if (vertIndx[0] != newVertCount
            && vertIndx[1] != newVertCount
            && vertIndx[2] != newVertCount) {
            for (unsigned int j = 0; j < 3; j++) {
                trianglesList.push_back(vertIndx[j]);
            }
            newTriangleCount++;
        }
    }
    unsigned int* triangles = new unsigned int[newTriangleCount * 3]{ 0 };
    std::copy(trianglesList.begin(), trianglesList.end(), triangles);

    // Redo the arrays, so that every triangle now has it own 3 vertices 
    if (altColAndMesh) {
        float* verticesAlt = new float[newTriangleCount * 9]{ 0 };
        float* normalsAlt = new float[newTriangleCount * 9]{ 0 };
        unsigned char* coloursAlt = new unsigned char[newTriangleCount * 9]{ 0 };

        unsigned int* atomIndexAlt = new unsigned int[newTriangleCount * 3]{ 0 };
        float* valuesAlt = new float[newTriangleCount * 3]{ 0 };
        unsigned int* trianglesAlt = new unsigned int[newTriangleCount * 3]{ 0 };

        for (unsigned int i = 0; i < newTriangleCount; i++) {   // i triangles
            for (uint8_t j = 0; j < 3; j++) {          // j vertices
                for (uint8_t k = 0; k < 3; k++) {      // k koordinates
                    verticesAlt[i * 9 + j * 3 + k] = vertices[triangles[i * 3 + j] * 3 + k];
                    normalsAlt[i * 9 + j * 3 + k] = normals[triangles[i * 3 + j] * 3 + k];
                    //coloursAlt[i * 9 + j * 3 + k] = colours[triangles[i * 3 + j] * 3 + k];    // the original colors
                    coloursAlt[i * 9 + j * 3 + k] = coloringFunction(i, k);
                }

                atomIndexAlt[i * 3 + j] = atomIndex[triangles[i * 3 + j]];
                valuesAlt[i * 3 + j] = values[triangles[i * 3 + j]];
                trianglesAlt[i * 3 + j] = i * 3 + j;
            }
        }

        // initialize the new mesh with the cut data
        cutMesh->SetVertexData(newTriangleCount * 3, verticesAlt, normalsAlt, coloursAlt, NULL, true);
        cutMesh->SetTriangleData(newTriangleCount, trianglesAlt, true);
        cutMesh->AddVertexAttribPointer(atomIndexAlt);
        cutMesh->AddVertexAttribPointer(valuesAlt);
        cutMesh->SetMaterial(NULL);

        delete[] vertices, normals, colours, atomIndex, values, triangles;
    }
    else {
        // initialize the new mesh with the cut data
        cutMesh->SetVertexData(newVertCount, vertices, normals, colours, NULL, true);
        cutMesh->SetTriangleData(newTriangleCount, triangles, true);
        cutMesh->AddVertexAttribPointer(atomIndex);
        cutMesh->AddVertexAttribPointer(values);
        cutMesh->SetMaterial(NULL);
    }

    delete[] oldVertIndices, triangles;

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

float* ViewOptimizationRenderer::removeDuplicatVertices(const float* vertices, const unsigned int arrLen, const unsigned int vertexCount) {
    glm::vec3 currentVert = glm::vec3(0);
    std::list<glm::vec3> vecList;

    // save the vertex array in alternate form
    for (unsigned int i = 0; i < arrLen; i++) {
        currentVert = glm::vec3(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2]);
        vecList.push_back(currentVert);
    }

    // remove duplicates
    vecList.unique();

    std::cout << "VEC LIST SIZE: " << vecList.size() << "\n";

    // turn the list of glm::vec3 back into an array of float
    float* uniqueVertices = new float[vecList.size() * 3]{ 0 };
    for (unsigned int i = 0; i < vecList.size(); i++) {
        uniqueVertices[i * 3] = vecList.front().x;
        uniqueVertices[i * 3 + 1] = vecList.front().y;
        uniqueVertices[i * 3 + 2] = vecList.front().z;
        vecList.pop_front();
    }

    return uniqueVertices;
}


char ViewOptimizationRenderer::coloringFunction(const unsigned int i, const uint8_t k) {
    switch (k) {
        case 0: // Red
            return char((i % 255) + 1);
            break;
        case 1: // Green
            return char(int(i / 255) % 255 + 1);
            break;
        case 2: //Blue
            return char((int(i / (255 * 255)) % 255 + 1));
            break;
        default:
            return char(0);
    };
}

unsigned int* ViewOptimizationRenderer::pixelColorCounter(const uint8_t* textureData, const unsigned int triangleCount, const unsigned int textureDataLength) {
    unsigned int* pixelForEachTriangle = new unsigned int[triangleCount + 1] {0};
    unsigned int triangleIndex;

    // works inverse of 'coloringFunction', decoding the rgb values to triangle indices
    for (unsigned int i = 0; i < textureDataLength / 3; i++) {
        if (textureData[i * 3] == 0 && textureData[i * 3 + 1] == 0 && textureData[i * 3 + 2] == 0) {
            pixelForEachTriangle[triangleCount]++;
        }
        else {
            triangleIndex = (textureData[i * 3] - 1) +
                (textureData[i * 3 + 1] - 1) * 255 +
                (textureData[i * 3 + 2] - 1) * 255 ^ 2;
            if (triangleIndex <= triangleCount) {
                pixelForEachTriangle[triangleIndex]++;
            }
        }
    }

    return pixelForEachTriangle;
}

float* ViewOptimizationRenderer::evaluateViewpoints(
    mmstd_gl::CallRender3DGL& call, const float* vertices, const unsigned int vertexCount, const glm::vec3 ligandCenter, const float radius) {

    /* ------- Setup variables ------- */

    float* directionVE = new float[vertexCount] {0};
    glm::vec3 direction = glm::vec3(0);
    std::shared_ptr<glowl::Texture2D> colorTexture = nullptr;
    megamol::compositing_gl::CallTexture2D* ct2d = NULL;

    unsigned int height = 0;
    unsigned int width = 0;
    unsigned int rgb = 3;
    unsigned int storageSize = 0;

    float a_t = 0;      // since the background also counts as a triangle, otherwise sum up visible triangel surface
    float a_i = 0;

    for (unsigned int i = 0; i < vertexCount; i++) {    // Iterate every vertex. Here, vertices are synonymous with view directions

        /* ------- reposition the camera ------- */

        direction = glm::vec3(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2]);
        direction = glm::normalize(direction);

        setCallCamera(call, direction, ligandCenter, radius, this->camDistFactor.Param<core::param::FloatParam>()->Value());

        /* ------- get the new texture ------- */

        ct2d = this->getTexture_.CallAs<megamol::compositing_gl::CallTexture2D>();
        if (ct2d != NULL) {
            (*ct2d)(0);
            colorTexture = ct2d->getData();
        }
        if (colorTexture == nullptr) {
            return directionVE;
        }

        height = colorTexture->getHeight();
        width = colorTexture->getWidth();
        storageSize = height * width * rgb;
        uint8_t* textureData = new uint8_t[storageSize];

        glGetTextureImage(colorTexture->getName(), 0, GL_RGB, GL_UNSIGNED_BYTE, storageSize, textureData);

        /* ------- calculate the viewpoint entropy ------- */

        unsigned int* pixelForEachTriangle = pixelColorCounter(textureData, currentTargetMeshData->GetTriCount(), storageSize);      // output will be + 1, due to the background


        a_t = height * width;      // since the background also counts as a triangle, this definition of a_t is used. Otherwise sum up all visible triangel surfaces

        if (a_t != 0) {
            for (unsigned int triangle = 0; triangle < currentTargetMeshData->GetTriCount() + 1; triangle++) {
                a_i = pixelForEachTriangle[triangle];

                if (a_i != 0) {
                    directionVE[i] += -a_i / a_t * log2(a_i / a_t);
                }
            }
        }
    }

    return directionVE;
}
