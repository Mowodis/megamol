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
    /*
    megamol::geocalls_gl::CallTriMeshDataGL* ctmd =
        this->getMacromoleculeMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();
    */

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
    /* ============= Simple: set camera position ============= */

    if (isInputXPosChanged) {
        // Calculate naive camera position based on ligand center position

        // get pointer to MolecularDataCall
        MolecularDataCall* mol = this->getLigandPDBData_.CallAs<MolecularDataCall>();
        if (mol == NULL)
            return false;

        // set call time
        //mol->SetCalltime(callTime);
        // set frame ID and call data
        //mol->SetFrameID(static_cast<int>(callTime));

        if (!(*mol)(protein_calls::MolecularDataCall::CallForGetData))
            return false;
        // check if atom count is zero
        if (mol->AtomCount() == 0)
            return true;

        float* pos0 = new float[mol->AtomCount() * 3];

        memcpy(pos0, mol->AtomPositions(), mol->AtomCount() * 3 * sizeof(float));

        // get ligand ceneter coordinates
        glm::vec3 ligandCenter = glm::vec3(0, 0, 0);
        for (unsigned int i = 0; i < mol->AtomCount() * 3; i++) {
            ligandCenter[i % 3] += pos0[i] / mol->AtomCount();
        }

        float bbEdgeLengts[] = {0.3, 0.4, 0.5};
        //cout << meanCord[1] << meanCord[2] << meanCord[3];
        

        // Change the camera position coordinates of the CallRenderer3DGL call
        auto& cam = call.GetCamera();
        Camera::Pose newCamPose = cam.getPose();  // have the old camera as default

        glm::vec3 camShift = glm::vec3(30,30,30);

        //newCamPose.position = glm::vec3(27, 16, 120);
        newCamPose.position = ligandCenter + camShift;
        newCamPose.direction = -camShift;

        const auto& cam_intrinsics = cam.get<Camera::PerspectiveParameters>();
        

        call.SetCamera(Camera(newCamPose, cam_intrinsics));



        isInputXPosChanged = false;

        delete[] pos0;
    }

    return true;
}

int temoFunction() {
    return 5;
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
