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

#include "compositing_gl/CompositingCalls.h"
//#include "compositing_gl/SimpleRenderTarget.h"     // Not located the in "Public Header Files"

using namespace megamol;
using namespace megamol::protein_gl;

//using namespace megamol::core::param;
//using namespace megamol::core::view;
//using namespace megamol::compositing_gl;


using namespace megamol::core::view;

using namespace std;

ViewOptimizationRenderer::ViewOptimizationRenderer()
        : getMacromoleculeMeshData_("getMacromoleculeMeshData",
              "Connects the renderer to a data provider to retrieve macromolecule mesh data")
        , getLigandPDBData_("getLigandPDBData", "Connects the renderer to a data provider to retrieve lignad mesh data")
        , camXCoord_("camXCoordinate", "changes the x-coordinate of the camera positon")
        , camera_("Camera", "Access the latest camera snapshot")
        , exeCounter_(0)
        , version_(0) {

    getMacromoleculeMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    getMacromoleculeMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_OPTIONAL);
    this->MakeSlotAvailable(&getMacromoleculeMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
    this->MakeSlotAvailable(&this->getLigandPDBData_);

    this->camXCoord_.SetParameter(new megamol::core::param::FloatParam(0.0f));
    this->MakeSlotAvailable(&this->camXCoord_);

    //this->camera_.SetCallback(megamol::compositing_gl::CallCamera::ClassName(), "GetData", &ViewOptimizationRenderer::getCameraSnapshot);
    //this->camera_.SetCallback(megamol::compositing_gl::CallCamera::ClassName(), "GetMetaData", &ViewOptimizationRenderer::getMetaDataCallback);
    //this->MakeSlotAvailable(&this->camera_);
}

ViewOptimizationRenderer::~ViewOptimizationRenderer() {
    this->Release();
}

bool ViewOptimizationRenderer::create() {
    return true;
}

void ViewOptimizationRenderer::release() {}

bool ViewOptimizationRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    ++version_;
    last_used_camera_ = call.GetCamera();

    megamol::geocalls_gl::CallTriMeshDataGL* ctmd =
        this->getMacromoleculeMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();

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

bool ViewOptimizationRenderer::getCameraSnapshot(core::Call& caller) {
    auto cc = dynamic_cast<megamol::compositing_gl::CallCamera*>(&caller);

    if (cc == NULL)
        return false;

    cc->setData(last_used_camera_, version_);

    return true;
}

bool ViewOptimizationRenderer::getMetaDataCallback(core::Call& caller) {
    return true;
}
