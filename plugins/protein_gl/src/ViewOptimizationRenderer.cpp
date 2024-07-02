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
        , optimizeCamera_("OptimizeCamera", "Acts as a button to set the camera to view the ligand binding site") {

    getMacromoleculeMeshData_.SetCompatibleCall<megamol::geocalls_gl::CallTriMeshDataGLDescription>();
    getMacromoleculeMeshData_.SetNecessity(megamol::core::AbstractCallSlotPresentation::Necessity::SLOT_OPTIONAL);
    this->MakeSlotAvailable(&getMacromoleculeMeshData_);

    this->getLigandPDBData_.SetCompatibleCall<protein_calls::MolecularDataCallDescription>();
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
    /* ============= Simple: set camera position ============= */

    if (this->optimizeCamera_.Param<core::param::BoolParam>()->Value()) {
        // Calculate naive camera position based on ligand center position

        /* ------- Get Ligand Data ------- */
        // get pointer to MolecularDataCall
        MolecularDataCall* ligand = this->getLigandPDBData_.CallAs<MolecularDataCall>();
        if (ligand == NULL)
            return false;

        if (!(*ligand)(protein_calls::MolecularDataCall::CallForGetData))
            return false;
        // check if atom count is zero
        if (ligand->AtomCount() == 0)
            return true;

        /* ------- Get Macromolecule Data ------- */
        geocalls_gl::CallTriMeshDataGL* macroMol =
            this->getMacromoleculeMeshData_.CallAs<megamol::geocalls_gl::CallTriMeshDataGL>();
        if (macroMol == nullptr) {
            return false;
        }
        macroMol->SetFrameID(static_cast<int>(call.Time()));
        if (!(*macroMol)(1)) {
            return false;
        }

        macroMol->SetFrameID(static_cast<int>(call.Time()));
        if (!(*macroMol)(0)) {
            return false;
        }

        /* ------- Make Camera Calculations ------- */

        // get ligand ceneter coordinates
        float* pos0 = new float[ligand->AtomCount() * 3];

        memcpy(pos0, ligand->AtomPositions(), ligand->AtomCount() * 3 * sizeof(float));

        glm::vec3 ligandCenter = glm::vec3(0, 0, 0);
        for (unsigned int i = 0; i < ligand->AtomCount() * 3; i++) {
            ligandCenter[i % 3] += pos0[i] / ligand->AtomCount();
        }

        // get radius, which will be used to select vertices around the ligand center
        // i.e get max atom distance from ligand center + atom sphere radius + gap of macro molecule-ligand
        float radius = 0;
        for (int i = 0; i < ligand->AtomCount(); i++) {
            float dist = glm::distance(glm::vec3(pos0[i * 3], pos0[i * 3 + 1], pos0[i * 3 + 2]), ligandCenter);
            if (dist > radius) {
                radius = dist;
            }
        }
        // adding value approximatly of average gap between ligand and macro molecule (arbitrarily chosen)
        // changes in this value might/should/will significantly change the camer orientation
        radius += 2;

        // get obj: the (first) object of macromolecule caller slot. Used for the vertex data
        const auto& obj = macroMol->Objects()[0];
        auto vertCount = obj.GetVertexCount();
        auto triCount = obj.GetTriCount();

        // find all macro molecule tie mesh vertices that are within the 'radius' aroud the ligand center
        // std::list<unsigned int> selectedVerticesIndex;
        glm::vec3 naiveCamDirection = glm::vec3(0, 0, 0);
        unsigned int nrSelectedVertices = 0;
        for (int i = 0; i < vertCount; i++) {
            glm::vec3 currentVertex = glm::vec3(obj.GetVertexPointerFloat()[i * 3],
                obj.GetVertexPointerFloat()[i * 3 + 1],
                obj.GetVertexPointerFloat()[i * 3 + 2]);
            if (glm::distance(currentVertex, ligandCenter) <= radius) {
                naiveCamDirection += glm::normalize(ligandCenter - currentVertex);
                nrSelectedVertices++;
            }
        }
        naiveCamDirection = glm::normalize(naiveCamDirection);

        // info and test outputs
        cout << "Number of macro molecule object: " << macroMol->Count() << "\n";
        cout << "Vertex count: " << vertCount << "\n";
        cout << "Tri(angle) count: " << triCount << "\n";
        cout << "Vertex selection radius: " << radius << "\n";
        cout << "Nr selected vertices: " << nrSelectedVertices << "\n";
        cout << "Average direction: " << naiveCamDirection.x << ", " << naiveCamDirection.y << ", " << naiveCamDirection.z << "\n";
        cout << "Normalization test: " << glm::normalize(glm::vec3(-4, 1, 0)).x << ", " << glm::normalize(glm::vec3(-4, 1, 0)).y << ", " << glm::normalize(glm::vec3(-4, 1, 0)).z << " \n ";
        cout << "Vector Length test: " << glm::vec3(4, 1, 0).length() << "\n";
        //cout << "Color test: " << obj.GetColourPointerFloat()[0] << "\n";

        /* ------- Set New Camera ------- */

        // Change the camera position coordinates of the CallRenderer3DGL call
        auto& cam = call.GetCamera();
        Camera::Pose newCamPose = cam.getPose();  // have the old camera as default

        // camera-ligandcenter distance is _ times the vertex selection radius
        float camDistFactor = 16.0f;
        newCamPose.position = ligandCenter + naiveCamDirection * glm::vec3(radius * camDistFactor);
        newCamPose.direction = -naiveCamDirection;

        const auto& cam_intrinsics = cam.get<Camera::PerspectiveParameters>();
        
        call.SetCamera(Camera(newCamPose, cam_intrinsics));

        this->optimizeCamera_.Param<core::param::BoolParam>()->SetValue(false);
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
