/**
 * MegaMol
 * Copyright (c) 2024, MegaMol Dev Team
 * Author: Marcel Heine (Bachelor Student)
 * All rights reserved.
 */

#pragma once

#include "mmcore/CalleeSlot.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"
#include "mmstd_gl/renderer/Renderer3DModuleGL.h"

namespace megamol::mmstd_gl {

/*
* CameraOverwrite class
*/
class CameraOverwrite : public megamol::mmstd_gl::Renderer3DModuleGL {

public:

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "CameraOverwrite";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Overwrites the 'Camera' parameter of a 'CallRender3DGL' call.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    /** Ctor. */
    CameraOverwrite();

    /** Dtor. */
    ~CameraOverwrite() override;

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    /**
     * Implementation of 'release'.
     */
    void release() override;

private:
    /**
     * The get extents callback. The module should set the members of
     * 'call' to tell the caller the extents of its data (bounding boxes
     * and times).
     *
     * @param call The calling call.
     *
     * @return The return value of the function.    
     */
    bool GetExtents(mmstd_gl::CallRender3DGL& call) override;

    /**
     * The Open GL Render callback.
     *
     * @param call The calling call.
     * @return The return value of the function.
     */
    bool Render(mmstd_gl::CallRender3DGL& call) override;

    /**
     * Camera input slot with which the 'call' camera is being replaced with 
     */
    core::CallerSlot getCamera_;
};

}
