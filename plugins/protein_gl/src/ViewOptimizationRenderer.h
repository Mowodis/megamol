/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team
 * Author: Marcel Heine (Bachelor Student)
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

    /*
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "ViewOptimizationRenderer";
    }

    /*
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Camera viwepoint optimizer for protein-ligand docking scenes";
    }
    /*
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    /* Ctor. */
    ViewOptimizationRenderer();

    /* Dtor. */
    ~ViewOptimizationRenderer() override;

protected:
    /*
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    /*
     * Implementation of 'release'.
     */
    void release() override;

    /*
      * The get extents callback. The module should set the members of
      * 'call' to tell the caller the extents of its data (bounding boxes
      * and times).
      *
      * @param call The calling call.
      *
      * @return The return value of the function.
      */
    bool GetExtents(mmstd_gl::CallRender3DGL& call) override;

    /*
     * The Open GL Render callback.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    bool Render(mmstd_gl::CallRender3DGL& call) override;

    /* =========== Global Variables =========== */

    /* The bounding box */
    vislib::math::Cuboid<float> bbox;

    /* Color relode variable */
    bool reload_colors;

    /* The data update hash */
    SIZE_T datahash;

    /* The camera position currently used during viewpoint sampling via viewpoint entropy*/
    mmstd_gl::CallRender3DGL currentSamplingCamera;

private:
    /* Require the mesh data from protein and ligand, as well as a texture*/
    core::CallerSlot getTargetMeshData_;
    core::CallerSlot getLigandPDBData_;
    core::CallerSlot getTexture_;

    /* Output of the targets cut triangel mesh data */
    core::CalleeSlot _cutTriangleMesh;

    /* Output a 'CallRender3DGL' object, containing only a camera */
    core::CalleeSlot _sampleSphereCamera;


    /* Module parameters */
    core::param::ParamSlot optimizeCamera;
    core::param::ParamSlot refreshTargetMesh;
    /** parameter slot for mesh rendering mode */
    megamol::core::param::ParamSlot renderTargetMeshModeParam;
    /** MSMS detail parameter */
    megamol::core::param::ParamSlot molRadiusSummand;
    /** Distance of camera to ligand center when optimizing*/
    megamol::core::param::ParamSlot camDistFactor;

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

    /*
     * Relays a 'CallRender3DGL' object containung only the most recent camera for viewpoint entriopy sampling
     */
    bool getSamplingCamera(core::Call& caller);

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
     *
     * @param call : rendering call
     * @param direction : normalized vector
     * @param ligCenter : coordinates of the averaged atom positions of the ligand molecule
     * @param radius : furthes atom of the ligand + 'molRadiusSummand'
     * @camDistFactor : multiplyed with 'radius' to determine how far from the 'ligCenter' the camera is being placed 
     */
    void setCallCamera(megamol::mmstd_gl::CallRender3DGL& call, glm::vec3 direction, glm::vec3 ligCenter, float radius, float camDistFactor);

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
     * Removes duplicate vertices from flaot array.
     * 
     *
     * @param vertices : float array, with sets of three coresponding to x,y,z coordinates
     * @param arrLen : length of 'vertices' array = 0 (mod 3)
     * @param vertexCount : number of unique vertices or sets of three
     */
    float* removeDuplicatVertices(const float* vertices, const unsigned int arrLen, const unsigned int vertexCount);

    /*
     * Returns a char values of 0 to 255 to be interpreted as rgb values.
     * Intended use seen in 'naiveCavetyCutter' as a way to implicitly encode the target mesh faces  
     *
     * @param i : 0 <= i <= 16581375, otherwise the colering is not unique anymore and loops under assumption that all posibilities for k have been exhaused as well.
     * @param k : 0 <= k <= 2, otherwise 0 is returned 
     */
    char coloringFunction(const unsigned int i, const uint8_t k);

    /*
     * Counts the occurences of uniquely colored pixels
     * Maps the pixel color to an indice the same way 'naiveCavetyCutter' maps triangle indices to colors when 'altColAndMesh' = true
     * Background color assumed black rgb = (0,0,0), stored in the last index.  
     *
     * @param textureData : rgb texture data, recolored by a cavetyCutter function
     * @param nrColors : number of colors or triangles + 1
     * @param´textureDataLength : index count aka lendth if the 'textureData' array 
     */
    unsigned int* pixelColorCounter(const uint8_t* textureData, const unsigned int triangleCount, const unsigned int textureDataLength);

    /*
     * Calculates the viewpoint entropy for every vertex at a set distance from the ligand center
     *
     * @param vertices : float array containing sets of three, interpreted as a vector with the ligand center, forming camera perspectives
     * @param vertexCount : length of 'vertices' divided by 3
     * @param distance : distance of the camera to the ligand center
     * @param ligandCenter : averaged atom positions of the ligand molecule
     */
    float* evaluateViewpoints(mmstd_gl::CallRender3DGL& call, const float* vertices, const unsigned int vertexCount, const glm::vec3 ligandCenter, const float radius);

    /*
     * Calculates the viewpoint entropy of the current scene using an icospheres vertices as sampling viewpoints
     */
    float calculateViewpointEntropy();
}; 

} 
