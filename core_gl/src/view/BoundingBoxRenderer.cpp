/*
 * BoundingBoxRenderer.cpp
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */
#include "stdafx.h"
#include "mmcore_gl/view/BoundingBoxRenderer.h"

#include "mmcore/CoreInstance.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/IntParam.h"

#include "vislib_gl/graphics/gl/ShaderSource.h"
#include "mmcore/utility/log/Log.h"

#include "glm/glm.hpp"
#include <glm/ext.hpp>
#include "glm/gtc/matrix_transform.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "mmcore_gl/utility/ShaderSourceFactory.h"

using namespace megamol::core_gl;
using namespace megamol::core_gl::view;

/*
 * BoundingBoxRenderer::BoundingBoxRenderer
 */
BoundingBoxRenderer::BoundingBoxRenderer(void)
    : RendererModule<CallRender3DGL, ModuleGL>()
    , enableBoundingBoxSlot("enableBoundingBox", "Enables the rendering of the bounding box")
    , boundingBoxColorSlot("boundingBoxColor", "Color of the bounding box")
    , smoothLineSlot("smoothLines", "Enables the smoothing of lines (may look strange on some setups)")
    , enableViewCubeSlot("enableViewCube", "Enables the rendering of the view cube")
    , viewCubePosSlot("viewCubePosition", "Position of the view cube")
    , viewCubeSizeSlot("viewCubeSize", "Size of the view cube")
    , vbo(0)
    , ibo(0)
    , va(0)
    , boundingBoxes () {

    this->enableBoundingBoxSlot.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->enableBoundingBoxSlot);

    this->boundingBoxColorSlot.SetParameter(new core::param::ColorParam("#ffffffff"));
    this->MakeSlotAvailable(&this->boundingBoxColorSlot);

    this->smoothLineSlot.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->smoothLineSlot);

    this->enableViewCubeSlot.SetParameter(new core::param::BoolParam(false));
    this->MakeSlotAvailable(&this->enableViewCubeSlot);

    this->viewCubePosSlot.SetParameter(new core::param::EnumParam(3));
    this->viewCubePosSlot.Param<core::param::EnumParam>()->SetTypePair(0, "bottom left");
    this->viewCubePosSlot.Param<core::param::EnumParam>()->SetTypePair(1, "bottom right");
    this->viewCubePosSlot.Param<core::param::EnumParam>()->SetTypePair(2, "top left");
    this->viewCubePosSlot.Param<core::param::EnumParam>()->SetTypePair(3, "top right");
    this->MakeSlotAvailable(&this->viewCubePosSlot);

    this->viewCubeSizeSlot.SetParameter(new core::param::IntParam(100));
    this->MakeSlotAvailable(&this->viewCubeSizeSlot);

    this->MakeSlotAvailable(&this->chainRenderSlot);
    this->MakeSlotAvailable(&this->renderSlot);
}

/*
 * BoundingBoxRenderer::~BoundingBoxRenderer
 */
BoundingBoxRenderer::~BoundingBoxRenderer(void) { this->Release(); }

/*
 * BoundingBoxRenderer::create
 */
bool BoundingBoxRenderer::create(void) {
    // TODO the vislib shaders have to die a slow and painful death
    vislib_gl::graphics::gl::ShaderSource bbVertSrc, bbFragSrc;
    auto ssf = std::make_shared<core_gl::utility::ShaderSourceFactory>(instance()->Configuration().ShaderDirectories());
    if (!ssf->MakeShaderSource("boundingbox::vertex", bbVertSrc)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to load vertex shader source for bounding box line shader");
    }
    if (!ssf->MakeShaderSource("boundingbox::fragment", bbFragSrc)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to load fragment shader source for bounding box line shader");
    }
    try {
        if (!this->lineShader.Create(bbVertSrc.Code(), bbVertSrc.Count(), bbFragSrc.Code(), bbFragSrc.Count())) {
            throw vislib::Exception("Shader creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to create bounding box line shader: %s\n", e.GetMsgA());
        return false;
    }

    vislib_gl::graphics::gl::ShaderSource vcVertSrc, vcFragSrc;
    if (!ssf->MakeShaderSource("viewcube::vertex", vcVertSrc)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to load vertex shader source for view cube shader");
    }
    if (!ssf->MakeShaderSource("viewcube::fragment", vcFragSrc)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to load fragment shader source for view cube shader");
    }
    try {
        if (!this->cubeShader.Create(vcVertSrc.Code(), vcVertSrc.Count(), vcFragSrc.Code(), vcFragSrc.Count())) {
            throw vislib::Exception("Shader creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("Unable to create view cube shader: %s\n", e.GetMsgA());
        return false;
    }

    // the used vertex order resembles the one used in the old View3D
    std::vector<float> vertexCoords = {-1.0f, -1.0f, -1.0f, 1.0f, -1.0f, -1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 1.0f, -1.0f,
        -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    std::vector<unsigned int> vertexIndices = {0, 2, 3, 1, 4, 5, 7, 6, 2, 6, 7, 3, 0, 1, 5, 4, 0, 4, 6, 2, 1, 3, 7, 5};

    glCreateBuffers(1, &this->vbo);
    glCreateBuffers(1, &this->ibo);
    glCreateVertexArrays(1, &this->va);

    glBindVertexArray(this->va);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->ibo);

    glEnableVertexAttribArray(0);

    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertexCoords.size(), vertexCoords.data(), GL_STATIC_DRAW);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * vertexIndices.size(), vertexIndices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    return true;
}

/*
 * BoundingBoxRenderer::release
 */
void BoundingBoxRenderer::release(void) {
    this->lineShader.Release();
    if (this->va != 0) {
        glDeleteVertexArrays(1, &this->va);
        this->va = 0;
    }
    if (this->vbo != 0) {
        glDeleteBuffers(1, &this->vbo);
        this->vbo = 0;
    }
    if (this->ibo != 0) {
        glDeleteBuffers(1, &this->ibo);
        this->ibo = 0;
    }
}

/*
 * BoundingBoxRenderer::GetExtents
 */
bool BoundingBoxRenderer::GetExtents(CallRender3DGL& call) {

    CallRender3DGL* chainedCall = this->chainRenderSlot.CallAs<CallRender3DGL>();
    if (chainedCall != nullptr) {
        *chainedCall = call;
        if ((*chainedCall)(core::view::AbstractCallRender::FnGetExtents)) {
            call = *chainedCall;
            this->boundingBoxes = call.AccessBoundingBoxes();
            return true;
        }
    }
    megamol::core::utility::log::Log::DefaultLog.WriteError(
        "The BoundingBoxRenderer does not work without a renderer attached to its right");
    return false;
}

/*
 * BoundingBoxRenderer::Render
 */
bool BoundingBoxRenderer::Render(CallRender3DGL& call) {

    core::view::Camera cam = call.GetCamera();
    auto const lhsFBO = call.GetFramebuffer();

    glm::mat4 view = cam.getViewMatrix();
    glm::mat4 proj = cam.getProjectionMatrix();
    glm::mat4 mvp = proj * view;

    CallRender3DGL* chainedCall = this->chainRenderSlot.CallAs<CallRender3DGL>();
    if (chainedCall == nullptr) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "The BoundingBoxRenderer does not work without a renderer attached to its right");
        return false;
    }

    auto smoothLines = this->smoothLineSlot.Param<core::param::BoolParam>()->Value();

    lhsFBO->bind();
    glViewport(0, 0, lhsFBO->getWidth(), lhsFBO->getHeight());

    bool renderRes = true;
    if (this->enableBoundingBoxSlot.Param<core::param::BoolParam>()->Value()) {
        renderRes &= this->RenderBoundingBoxBack(mvp, this->boundingBoxes, smoothLines);
    }

    //glBindFramebuffer(GL_FRAMEBUFFER, 0);

    *chainedCall = call;
    renderRes &= (*chainedCall)(core::view::AbstractCallRender::FnRender);

    lhsFBO->bind();
    glViewport(0, 0, lhsFBO->getWidth(), lhsFBO->getHeight());

    if (this->enableBoundingBoxSlot.Param<core::param::BoolParam>()->Value()) {
        renderRes &= this->RenderBoundingBoxFront(mvp, this->boundingBoxes, smoothLines);
    }
    if (this->enableViewCubeSlot.Param<core::param::BoolParam>()->Value()) {
        renderRes &= this->RenderViewCube(call);
    }

    // glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return renderRes;
}

/*
 * BoundingBoxRenderer::RenderBoundingBoxFront
 */
bool BoundingBoxRenderer::RenderBoundingBoxFront(const glm::mat4& mvp, const core::BoundingBoxes_2& bb, bool smoothLines) {
    glm::vec3 bbmin = glm::vec3(bb.BoundingBox().Left(), bb.BoundingBox().Bottom(), bb.BoundingBox().Back());
    glm::vec3 bbmax = glm::vec3(bb.BoundingBox().Right(), bb.BoundingBox().Top(), bb.BoundingBox().Front());

    auto colptr = this->boundingBoxColorSlot.Param<core::param::ColorParam>()->Value();

    this->lineShader.Enable();

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    if (smoothLines) glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.75f);
    glPolygonMode(GL_FRONT, GL_LINE);

    glBindVertexArray(this->va);
    glEnableVertexAttribArray(0);

    glUniformMatrix4fv(this->lineShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
    glUniform3fv(this->lineShader.ParameterLocation("bbMin"), 1, glm::value_ptr(bbmin));
    glUniform3fv(this->lineShader.ParameterLocation("bbMax"), 1, glm::value_ptr(bbmax));
    glUniform3f(this->lineShader.ParameterLocation("color"), colptr[0], colptr[1], colptr[2]);

    glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);

    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    if (smoothLines) glDisable(GL_LINE_SMOOTH);
    glLineWidth(1.0f);
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT, GL_FILL);

    this->lineShader.Disable();

    return true;
}

/*
 * BoundingBoxRenderer::RenderBoundingBoxBack
 */
bool BoundingBoxRenderer::RenderBoundingBoxBack(const glm::mat4& mvp, const core::BoundingBoxes_2& bb, bool smoothLines) {
    glm::vec3 bbmin = glm::vec3(bb.BoundingBox().Left(), bb.BoundingBox().Bottom(), bb.BoundingBox().Back());
    glm::vec3 bbmax = glm::vec3(bb.BoundingBox().Right(), bb.BoundingBox().Top(), bb.BoundingBox().Front());

    auto color = this->boundingBoxColorSlot.Param<core::param::ColorParam>()->Value();

    this->lineShader.Enable();

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    if (smoothLines) glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.25f);
    glPolygonMode(GL_BACK, GL_LINE);

    glBindVertexArray(this->va);
    glEnableVertexAttribArray(0);

    glUniformMatrix4fv(this->lineShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
    glUniform3fv(this->lineShader.ParameterLocation("bbMin"), 1, glm::value_ptr(bbmin));
    glUniform3fv(this->lineShader.ParameterLocation("bbMax"), 1, glm::value_ptr(bbmax));
    glUniform3f(this->lineShader.ParameterLocation("color"), color[0], color[1], color[2]);

    glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);

    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    if (smoothLines) glDisable(GL_LINE_SMOOTH);
    glLineWidth(1.0f);
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_BACK, GL_FILL);

    this->lineShader.Disable();

    return true;
}

/*
 * BoundingBoxRenderer::RenderViewCube
 */
bool BoundingBoxRenderer::RenderViewCube(CallRender3DGL& call) {
    // Get camera orientation
    auto cam = call.GetCamera();

    const auto rotation = glm::mat4(glm::inverse(glm::mat3(cam.getViewMatrix())));

    // Create view/model and projection matrices
    const float dist = 2.0f / std::tan(glm::radians(30.0f) / 2.0f);

    glm::mat4 model(1.0f);
    model[3][2] = -dist;

    const auto proj = glm::perspective(glm::radians(30.0f), 1.0f, 0.1f, 100.0f);

    // Set state
    const auto depth_test = glIsEnabled(GL_DEPTH_TEST);
    if (depth_test) glDisable(GL_DEPTH_TEST);

    const auto culling = glIsEnabled(GL_CULL_FACE);
    if (!culling) glEnable(GL_CULL_FACE);

    // Set viewport
    std::array<GLint, 4> viewport;
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    const auto position = this->viewCubePosSlot.Param<core::param::EnumParam>()->Value();
    const auto size = this->viewCubeSizeSlot.Param<core::param::IntParam>()->Value();

    int x, y;

    switch (position) {
    case 0:
        x = viewport[0];
        y = viewport[1];
        break;
    case 1:
        x = viewport[2] - size;
        y = viewport[1];
        break;
    case 2:
        x = viewport[0];
        y = viewport[3] - size;
        break;
    case 3:
    default:
        x = viewport[2] - size;
        y = viewport[3] - size;
        break;
    }

    glViewport(x, y, size, size);

    // Render view cube
    this->cubeShader.Enable();

    glUniformMatrix4fv(this->cubeShader.ParameterLocation("rot_mx"), 1, false, glm::value_ptr(rotation));
    glUniformMatrix4fv(this->cubeShader.ParameterLocation("model_mx"), 1, false, glm::value_ptr(model));
    glUniformMatrix4fv(this->cubeShader.ParameterLocation("proj_mx"), 1, false, glm::value_ptr(proj));

    glDrawArrays(GL_TRIANGLES, 0, 36);

    this->cubeShader.Disable();

    // Restore viewport
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

    // Restore state
    if (depth_test) glEnable(GL_DEPTH_TEST);
    if (!culling) glDisable(GL_CULL_FACE);

    return true;
}
