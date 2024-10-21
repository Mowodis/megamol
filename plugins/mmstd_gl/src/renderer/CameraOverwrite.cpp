/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team
 * Author: Marcel Heine (Bachelor Student)
 * All rights reserved.
 */

#include "mmstd_gl/renderer/CameraOverwrite.h"

using namespace megamol::mmstd_gl;

CameraOverwrite::CameraOverwrite()
        : mmstd_gl::Renderer3DModuleGL()
        , getCamera_("getCamera", "Camera input slot, with which the camera in the 'chainRendering' call is being replaced with") {
    this->getCamera_.SetCompatibleCall<megamol::mmstd_gl::CallRender3DGLDescription>();
    this->getCamera_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&getCamera_);
}

CameraOverwrite::~CameraOverwrite() {
    this->Release();
}

bool CameraOverwrite::create() {
    return true;
}

void CameraOverwrite::release() {}

bool CameraOverwrite::GetExtents(CallRender3DGL& call) {
    return true;
}

bool CameraOverwrite::Render(CallRender3DGL& call) {

    CallRender3DGL* cr = this->getCamera_.CallAs<CallRender3DGL>();

    if (cr == nullptr) {
        return false;
    }

    if (!(*cr)(0)) {
        //std::cout << "CR IS NULL!!! \n";
        return false;
    }

    /*
    if (!(*cr)(0)) {
        return false;
    }
    */

    megamol::core::view::Camera cam = cr->GetCamera();

    /*
    std::cout << "Camera Overwrite: "
              << "x:" << cam.getPose().direction.x << ", y:" << cam.getPose().direction.y
              << ", z:" << cam.getPose().direction.z;
    */
    /*
     if (call.GetCamera().getPose().direction.x != cam.getPose().direction.x,
        call.GetCamera().getPose().direction.y != cam.getPose().direction.y,
        call.GetCamera().getPose().direction.z != cam.getPose().direction.z) {

        std::cout << "Camera Overwrite: "
                  << "x:" << cam.getPose().direction.x << ", y:" << cam.getPose().direction.y
                  << ", z:" << cam.getPose().direction.z;
    }
    */
    

    call.SetCamera(cam);

    
    

    return true;
}

