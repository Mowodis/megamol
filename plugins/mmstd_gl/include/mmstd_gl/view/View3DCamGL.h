/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team, bachelor student Marcel Heine
 * All rights reserved.
 */

/**
 * Module "View3DGL" modefied with new caler slot to accept a camera
*/

#pragma once

#include <glowl/FramebufferObject.hpp>

#include "mmcore/view/CameraControllers.h"
#include "mmstd/view/BaseView.h"
#include "mmstd_gl/renderer/CallRenderViewGL.h"
#include "mmstd_gl/view/AbstractViewGL.h"

namespace megamol::mmstd_gl::view {

class View3DCamGL : public core::view::BaseView<CallRenderViewGL, core::view::Camera3DController, AbstractViewGL> {

public:
    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "View3DCamGL";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "New and improved 3D View Module";
    }

    /** Ctor. */
    View3DCamGL();

    /** Dtor. */
    ~View3DCamGL() override;

    ImageWrapper Render(double time, double instanceTime) override;

    ImageWrapper GetRenderingResult() const override;

    /**
     * Resizes the framebuffer object and calls base class function that sets camera aspect ratio if applicable.
     *
     * @param width The new width.
     * @param height The new height.
     */
    void Resize(unsigned int width, unsigned int height) override;

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

private:
    core::CallerSlot _cameraSlot;
};

} // namespace megamol::mmstd_gl::view
