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

using namespace std;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getTargetMeshData_("getTargetMeshData",
              "Connects the renderer to a data provider to retrieve target mesh data")
        , getLigandPDBData_("getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , optimizeCamera_("OptimizeCamera", "Acts as a button to set the camera to view the ligand binding site") {

    this->getTargetMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    this->getTargetMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&getTargetMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
    this->getLigandPDBData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->getLigandPDBData_);

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

        /* ------- Make Naive Camera Calculations ------- */

        // get ligand ceneter coordinates
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

        // get obj: the (first) object of target caller slot. Used for the vertex data
        const auto& obj = target->Objects()[0];
        auto vertCount = obj.GetVertexCount();
        auto triCount = obj.GetTriCount();

        // find all target tie mesh vertices that are within the 'radius' aroud the ligand center
        // std::list<unsigned int> selectedVerticesIndex;
        glm::vec3 naiveCamDirection = glm::vec3(0, 0, 0);
        unsigned int nrSelectedVertices = 0;
        for (int i = 0; i < vertCount; i++) {
            glm::vec3 currentVertex = glm::vec3(obj.GetVertexPointerFloat()[i * 3], obj.GetVertexPointerFloat()[i * 3 + 1], obj.GetVertexPointerFloat()[i * 3 + 2]);

            if (glm::distance(currentVertex, ligandCenter) <= radius) {
                naiveCamDirection += glm::normalize(ligandCenter - currentVertex);
                nrSelectedVertices++;
            }
        }
        naiveCamDirection = glm::normalize(naiveCamDirection);
        cout << "CORDS: " << naiveCamDirection[0] << "," << naiveCamDirection[1] << "," << naiveCamDirection[2];

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

glm::vec3 ViewOptimizationRenderer::moleculeCenter(float* positions, unsigned int atomCount) {
    glm::vec3 ligandCenter = glm::vec3(0, 0, 0);
    for (unsigned int i = 0; i < atomCount * 3; i++) {
        ligandCenter[i % 3] += positions[i] / atomCount;
    }

    return ligandCenter;
}

void ViewOptimizationRenderer::newCallCamera(mmstd_gl::CallRender3DGL& call, glm::vec3 direction, glm::vec3 ligCenter, float radius, float camDistFactor) {
    // Change the camera position coordinates of the CallRenderer3DGL call
    auto& cam = call.GetCamera();
    Camera::Pose newCamPose = cam.getPose();  // old camera as default

    // camera-ligandcenter distance is _ times the vertex selection radius
    newCamPose.position = ligCenter + direction * glm::vec3(radius * camDistFactor);
    newCamPose.direction = -direction;

    const auto& camIntrinsics = cam.get<Camera::PerspectiveParameters>();

    call.SetCamera(Camera(newCamPose, camIntrinsics));
}

//geocalls_gl::CallTriMeshDataGL::Mesh naiveCavetyCutter(geocalls_gl::CallTriMeshDataGL::Mesh mesh, glm::vec3 ligCenter) {

//}
