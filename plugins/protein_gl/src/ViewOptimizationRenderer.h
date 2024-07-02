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

namespace megamol::protein_gl {

class ViewOptimizationRenderer : public mmstd_gl::Renderer3DModuleGL {
public:
    static const char* ClassName() {
        return "ViewpointOptimizer";
    }

    static const char* Description() {
        return "Camera viwepoint optimizer for protein-ligand docking";
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

private:
    /* Require the mesh data from protein and ligand*/
    core::CallerSlot getMacromoleculeMeshData_;
    core::CallerSlot getLigandPDBData_;

    /* Module parameters */
    core::param::ParamSlot optimizeCamera_;

    bool isInputXPosChanged;

}; // class Vie//wpointOptimizer

} // namespace megamol::protein_gl
