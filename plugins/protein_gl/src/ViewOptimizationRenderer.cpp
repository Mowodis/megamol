/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team
 * All rights reserved.
 */

#include "ViewOptimizationRenderer.h"

// #include "mmcore/view/Camera.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"

//#include "geometry_calls_gl/CallTriMeshDataGL.h"



using namespace megamol;
using namespace megamol::protein_gl;

ViewOptimizationRenderer::ViewOptimizationRenderer() {

}

ViewOptimizationRenderer::~ViewOptimizationRenderer() {
    this->Release();
}

bool ViewOptimizationRenderer::create() {
    return true;
}

void ViewOptimizationRenderer::release() {}

bool ViewOptimizationRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    return true;
}

bool ViewOptimizationRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    return true;
}

