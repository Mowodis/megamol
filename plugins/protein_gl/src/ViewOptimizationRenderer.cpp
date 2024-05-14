/**
 * MegaMol
 * Copyright (c) 2024 by Marcel Heine
 * All rights reserved.
 */

#include "ViewOptimizationRenderer.h"

#include <glm/glm.hpp>

#include "geometry_calls_gl/CallTriMeshDataGL.h"
#include "mmcore/view/Camera.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"
#include "mmcore/param/FloatParam.h"

#include "protein_calls/MolecularDataCall.h"

using namespace megamol;
using namespace megamol::core::param;
using namespace megamol::core::view;
using namespace megamol::protein_gl;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getMacromoleculeMeshData_("getMacromoleculeMeshData",
              "Connects the renderer to a data provider to retrieve macromolecule mesh data")
        , getLigandPDBData_(
              "getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , camXCoord_("camXCoordinate", "changes the x-coordinate of the camera positon") {

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
    auto& cam = call.GetCamera();
    const auto& pos = cam.getPose();
    auto new_pose = pos;

    float x_pos = this->camXCoord_.Param<core::param::FloatParam>()->Value();
    auto const center_ = call.AccessBoundingBoxes().BoundingBox().CalcCenter();

    const glm::vec3 cam_pos = pos.position;
    const glm::vec3 new_cam_pos = glm::vec3(x_pos, center_.GetY(), center_.GetZ());

    new_pose.position = new_cam_pos;
    const auto& cam_intrinsics = cam.get<Camera::PerspectiveParameters>();

    call.SetCamera(Camera(new_pose, cam_intrinsics));


    return true;
}

bool ViewOptimizationRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    return true;
}
