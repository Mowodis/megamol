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



using namespace megamol::core::view;

using namespace std;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getMacromoleculeMeshData_("getMacromoleculeMeshData",
              "Connects the renderer to a data provider to retrieve macromolecule mesh data")
        , getLigandPDBData_("getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , camXCoord_("camXCoordinate", "changes the x-coordinate of the camera positon")
        , exeCounter_(0)
        , isInputXPosChanged(true) {

    getMacromoleculeMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    getMacromoleculeMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_OPTIONAL);
    this->MakeSlotAvailable(&getMacromoleculeMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->getLigandPDBData_);

    this->camXCoord_.SetParameter(new megamol::core::param::FloatParam(0.0f));
    this->MakeSlotAvailable(&this->camXCoord_);
}

ViewOptimizationRenderer::~ViewOptimizationRenderer() {
    this->Release();
}

bool ViewOptimizationRenderer::create() {
    return true;
}

void ViewOptimizationRenderer::release() {}

bool ViewOptimizationRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    megamol::geocalls_gl::CallTriMeshDataGL* ctmd =
        this->getMacromoleculeMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();

    /* 
    int tickTarget = 50;
    if (exeCounter_ == tickTarget) {
        cout << "Here's some info: \n";

        ctmd->SetFrameID(static_cast<int>(call.Time()));
        if (!(*ctmd)(1)) {
            return false;
        }

        ctmd->SetFrameID(static_cast<int>(call.Time()));
        if (!(*ctmd)(0)) {
            return false;
        }

        if (ctmd->Count() > 0) {
            cout << ctmd->Objects()[0].GetVertexCount() << "\n ";
        }

        exeCounter_++;
    }
    if (exeCounter_ < tickTarget) {
        exeCounter_++;
    }
    */
    /* ============= Dummy set camera x position ============= */

    if (isInputXPosChanged) {
        auto& cam = call.GetCamera();
        auto cam_pose = cam.getPose();

        //float x_pos = this->camXCoord_.Param<core::param::FloatParam>()->Value();
        //auto center_ = call.AccessBoundingBoxes().BoundingBox().CalcCenter();

        //glm::vec3 cam_pos = cam.getPose().position;
        glm::vec3 new_cam_pos = glm::vec3(27, 16, 120);

        cam_pose.position = new_cam_pos;
        const auto& cam_intrinsics = cam.get<Camera::PerspectiveParameters>();

        call.SetCamera(Camera(cam_pose, cam_intrinsics));

        isInputXPosChanged = false;
    }

    return true;
}

bool ViewOptimizationRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    auto ctmd = this->getMacromoleculeMeshData_.CallAs<geocalls_gl::CallTriMeshDataGL>();
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
