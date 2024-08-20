/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team, bachelor student Marcel Heine
 * All rights reserved.
 */

#pragma once

#include <glm/glm.hpp>

#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/view/Camera.h"
#include "mmstd_gl/renderer/Renderer3DModuleGL.h"
#include "mmcore/CalleeSlot.h"
#include "geometry_calls_gl/CallTriMeshDataGL.h"

namespace megamol::protein_gl {

class ViewOptimizationRenderer : public mmstd_gl::Renderer3DModuleGL {
public:
    static const char* ClassName() {
        return "ViewOptimizationRenderer";
    }

    static const char* Description() {
        return "Camera viwepoint optimizer for protein-ligand docking scenes";
    }

    static bool IsAvailable() {
        return true;
    }

    ViewOptimizationRenderer();
    ~ViewOptimizationRenderer() override;

protected:
    bool create() override;
    void release() override;
    bool GetExtents(mmstd_gl::CallRender3DGL& call) override;
    bool Render(mmstd_gl::CallRender3DGL& call) override;

    /** The bounding box */
    vislib::math::Cuboid<float> bbox;

    bool reload_colors;

    /** The data update hash */
    SIZE_T datahash;

private:
    /* Require the mesh data from protein and ligand*/
    core::CallerSlot getTargetMeshData_;
    core::CallerSlot getLigandPDBData_;

    /* Output of the targets cut triangel mesh data */
    core::CalleeSlot _cutTriangleMesh;

    /* Module parameters */
    core::param::ParamSlot optimizeCamera_;

    /* Stores the current target molecule object data */
    const geocalls_gl::CallTriMeshDataGL::Mesh* currentTargetMeshData;

    /* =========== Functions =========== */

    /*
     * Relays the mesh data of a target molecule to a 'ModernTrisoupRendererGL'
     */
    bool getDataCallback(core::Call& caller);


    /*
     * Extend function for molecular mesh data relay puproses
     */
    bool getExtentCallback(core::Call& caller);

    /*
     * Determine the center coordinate (x1,x2,x3) by averaging over all atom positions
     * defaults to (0,0,0) if no atoms are contained within the call
     *
     * Input: float array of length = 0 (mod 3), where sets of three denote the (x,y,z) coordinates of an atom
     */
    glm::vec3 moleculeCenter(float* coordinates, unsigned int atomCount);

    /*
     * Calculate a camera direction vector in a naive manner
     * Summs up the vectors from within range target mesh vertices to the ligand center and normalizes the result
     */
    glm::vec3 naiveCameraDirection(megamol::geocalls_gl::CallTriMeshDataGL* obj, const glm::vec3 center, const float radius);

    /*
     * Replace the camera position and direction of a call with new values
     * SIDEEFFECTS
     */
    void newCallCamera(megamol::mmstd_gl::CallRender3DGL& call, glm::vec3 direction, glm::vec3 ligCenter, float radius, float camDistFactor);

    /* remove all mesh vertices, that do not lay within a certain radius around a center coordinate
     * SIDEEFFECTS 
     */
    geocalls_gl::CallTriMeshDataGL::Mesh* naiveCavetyCutter(
        megamol::geocalls_gl::CallTriMeshDataGL::Mesh mesh, glm::vec3 centrioid, float radius);

    /*
     * checks if an element is contained within an array and returns its Index
     * If not found, returns the length of the array
     */
    unsigned int inArray(unsigned int* arr, unsigned int element, unsigned int arrSize);

}; // class Vie//wpointOptimizer

} // namespace megamol::protein_gl
