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
#include "protein_calls/MolecularDataCall.h"

namespace megamol::protein_gl {

class ViewOptimizationRenderer : public mmstd_gl::Renderer3DModuleGL {
public:

    /* Encodes the different working states of the target molecule triangle mesh */
    enum MeshRenderMode {
        WHOLE_MESH = 0,
        NAIVE_CAVETY = 1
    };

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
    /* Require the mesh data from protein and ligand, as well as a texture*/
    core::CallerSlot getTargetMeshData_;
    core::CallerSlot getLigandPDBData_;
    core::CallerSlot getTexture_;

    /* Output of the targets cut triangel mesh data */
    core::CalleeSlot _cutTriangleMesh;

    /* Module parameters */
    core::param::ParamSlot optimizeCamera;
    core::param::ParamSlot refreshTargetMesh;
    /** parameter slot for mesh rendering mode */
    megamol::core::param::ParamSlot renderTargetMeshModeParam;
    /** MSMS detail parameter */
    megamol::core::param::ParamSlot molRadiusSummand;

    /* Stores the current target molecule object data */
    const geocalls_gl::CallTriMeshDataGL::Mesh* currentTargetMeshData;

    /* Stores the mode that determines in which way the targets triangle mesh is being rendered */
    MeshRenderMode currentTargetMeshRenderMode;

    /* =========== Functions =========== */

    /*
     * Relays the mesh data of a target molecule to a 'ModernTrisoupRendererGL'
     */
    bool getDataCallback(core::Call& caller);


    /*
     * Extend function for molecular mesh data relay puproses
     */
    bool getExtentCallback(core::Call& caller);

    /**
     * Update parameter slots.
     */
    void UpdateParameters();

    /*
     * Determine the center coordinate (x1,x2,x3) by averaging over all atom positions
     * defaults to (0,0,0) if no atoms are contained within the call
     *
     * Input: float array of length = 0 (mod 3), where sets of three denote the (x,y,z) coordinates of an atom
     */
    glm::vec3 moleculeCenter(float* coordinates, unsigned int atomCount);


    /*
     * Finds the atom, that is furthest from the center.
     * The 'buffer' is added to that distance.
     */
    float moleculeRadius(megamol::protein_calls::MolecularDataCall* ligand, float* atomPos, glm::vec3 ligandCenter, float buffer = 0);

    /*
     * Calculate a camera direction vector in a naive manner
     * Summs up the vectors from within range target mesh vertices to the ligand center and normalizes the result
     */
    glm::vec3 naiveCameraDirection(megamol::geocalls_gl::CallTriMeshDataGL* obj, const glm::vec3 center, const float radius);

    /*
     * Replace the camera position and direction of a call with new values
     * SIDEEFFECTS on the 'call'
     */
    void newCallCamera(megamol::mmstd_gl::CallRender3DGL& call, glm::vec3 direction, glm::vec3 ligCenter, float radius, float camDistFactor);

    /*
     * remove all mesh vertices, that do not lay within a certain radius around a center coordinate
     */
    geocalls_gl::CallTriMeshDataGL::Mesh* naiveCavetyCutter(
        megamol::geocalls_gl::CallTriMeshDataGL::Mesh mesh, glm::vec3 centrioid, float radius, bool altColAndMesh);

    /*
     * checks if an element is contained within an array and returns its Index
     * If not found, returns the length of the array
     */
    unsigned int inArray(unsigned int* arr, unsigned int element, unsigned int arrSize);

    /*
     * Returns a char values of 0 to 255 to be interpreted as rgb values.
     * Intended use seen in 'naiveCavetyCutter' as a way to implicitly encode the target mesh faces  
     *
     * @param i : 0 <= i <= 16581375, otherwise the colering is not unique anymore and loops under assumption that all posibilities for k have been exhaused as well.
     * @param k : 0 <= k <= 2, otherwise 0 is returned 
     */
    char coloringFunction(unsigned int i, uint8_t k);

    /*
     * Counts the occurences of uniquely colored pixels
     * Maps the pixel color to an indice the same way 'naiveCavetyCutter' maps triangle indices to colors when 'altColAndMesh' = true
     * Background color assumed black rgb = (0,0,0), stored in the last index.  
     *
     * @param textureData : rgb texture data, recolored by a cavetyCutter function
     * @param nrColors : number of colors or triangles + 1
     * @param´textureDataLength : index count aka lendth if the 'textureData' array 
     */
    unsigned int* pixelColCounter(uint8_t* textureData, unsigned int triangleCount, unsigned int textureDataLength);


}; // class Vie//wpointOptimizer

} // namespace megamol::protein_gl
