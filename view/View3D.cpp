/*
 * View3D.cpp
 *
 * Copyright (C) 2008 - 2009 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "View3D.h"
#ifdef _WIN32
#include <windows.h>
#endif /* _WIN32 */
#include <GL/gl.h>
#include "view/CallCursorInput.h"
#include "view/CallRenderView.h"
#include "view/ClusterCameraParamOverride.h"
#include "param/BoolParam.h"
#include "param/ButtonParam.h"
#include "param/FloatParam.h"
#include "param/StringParam.h"
#include "CallRender.h"
#include "CallRender3D.h"
#include "utility/ColourParser.h"
#include "vislib/CameraParamsStore.h"
#include "vislib/Exception.h"
#include "vislib/Log.h"
#include "vislib/mathfunctions.h"
#include "vislib/Point.h"
#include "vislib/String.h"
#include "vislib/StringSerialiser.h"
#include "vislib/sysfunctions.h"
#include "vislib/Trace.h"
#include "vislib/Vector.h"

using namespace megamol::core;


/*
 * view::View3D::View3D
 */
view::View3D::View3D(void) : view::AbstractView(), cam(), camParams(),
        camOverrides(), setViewport(true), cursor2d(), modkeys(), rotator1(),
        rotator2(), zoomer1(), zoomer2(),
        rendererSlot("rendering", "Connects the view to a Renderer"),
        cursorInputSlot("cursorInput", "Slot for incoming cursor input"),
        lightDir(0.5f, -1.0f, -1.0f), isCamLight(true), bboxs(),
        animPlay("animPlay", "Bool parameter to play/stop the animation"),
        animSpeed("animSpeed", "Float parameter of animation speed in time frames per second"),
        showBBox("showBBox", "Bool parameter to show/hide the bounding box"),
        softCursor("softCursor", "Bool flag to activate software cursor rendering"),
        backCol("backCol", "The views background colour"),
        cameraSettingsSlot("camsettings", "The stored camera settings"),
        storeCameraSettingsSlot("storecam", "Triggers the storage of the camera settings"),
        restoreCameraSettingsSlot("restorecam", "Triggers the restore of the camera settings"),
        resetViewSlot("resetView", "Triggers the reset of the view"),
        timeFrame(0.0f), animTimer(true), fpsCounter(10), fpsOutputTimer(0),
        firstImg(false), frozenValues(NULL) {
    this->MakeSlotAvailable(view::AbstractView::getRenderViewSlot());

    this->camParams = this->cam.Parameters();
    this->camOverrides = new ClusterCameraParamOverride(this->camParams);

    vislib::graphics::ImageSpaceType defWidth(
        static_cast<vislib::graphics::ImageSpaceType>(100));
    vislib::graphics::ImageSpaceType defHeight(
        static_cast<vislib::graphics::ImageSpaceType>(100));

    this->camParams->SetVirtualViewSize(defWidth, defHeight);
    this->camParams->SetTileRect(vislib::math::Rectangle<float>(0.0f, 0.0f,
        defWidth, defHeight));

    this->bkgndCol[0] = 0.0f;
    this->bkgndCol[1] = 0.0f;
    this->bkgndCol[2] = 0.125f;

    this->rendererSlot.SetCompatibleCall<CallRenderDescription>();
    this->rendererSlot.SetCompatibleCall<CallRender3DDescription>();
    this->MakeSlotAvailable(&this->rendererSlot);

    this->cursorInputSlot.SetCallback(view::CallCursorInput::ClassName(), 
        view::CallCursorInput::FunctionName(0), &View3D::onSetCursor2DButtonState);
    this->cursorInputSlot.SetCallback(view::CallCursorInput::ClassName(), 
        view::CallCursorInput::FunctionName(1), &View3D::onSetCursor2DPosition);
    this->cursorInputSlot.SetCallback(view::CallCursorInput::ClassName(), 
        view::CallCursorInput::FunctionName(2), &View3D::onSetInputModifier);
    this->cursorInputSlot.SetCallback(view::CallCursorInput::ClassName(), 
        view::CallCursorInput::FunctionName(3), &View3D::onResetView);
    this->MakeSlotAvailable(&this->cursorInputSlot);

    // empty bounding box will trigger initialisation
    this->bboxs.Clear();

    // simple animation time controlling (TODO: replace)
    this->animPlay << new param::BoolParam(false);
    this->MakeSlotAvailable(&this->animPlay);
    this->animSpeed << new param::FloatParam(4.0f, 0.01f, 100.0f);
    this->MakeSlotAvailable(&this->animSpeed);
    this->showBBox << new param::BoolParam(true);
    this->MakeSlotAvailable(&this->showBBox);
    this->softCursor << new param::BoolParam(false, false);
    this->MakeSlotAvailable(&this->softCursor);

    this->backCol << new param::StringParam(utility::ColourParser::ToString(
        this->bkgndCol[0], this->bkgndCol[1], this->bkgndCol[2]));
    this->MakeSlotAvailable(&this->backCol);

    this->cameraSettingsSlot << new param::StringParam("");
    this->MakeSlotAvailable(&this->cameraSettingsSlot);

    this->storeCameraSettingsSlot << new param::ButtonParam(
        vislib::sys::KeyCode::KEY_MOD_ALT | vislib::sys::KeyCode::KEY_MOD_SHIFT | 'C');
    this->storeCameraSettingsSlot.SetUpdateCallback(&View3D::onStoreCamera);
    this->MakeSlotAvailable(&this->storeCameraSettingsSlot);

    this->restoreCameraSettingsSlot << new param::ButtonParam(
        vislib::sys::KeyCode::KEY_MOD_ALT | 'c');
    this->restoreCameraSettingsSlot.SetUpdateCallback(&View3D::onRestoreCamera);
    this->MakeSlotAvailable(&this->restoreCameraSettingsSlot);

    this->resetViewSlot << new param::ButtonParam(vislib::sys::KeyCode::KEY_HOME);
    this->resetViewSlot.SetUpdateCallback(&View3D::onResetView);
    this->MakeSlotAvailable(&this->resetViewSlot);

    this->ResetView();
}


/*
 * view::View3D::~View3D
 */
view::View3D::~View3D(void) {
    this->Release();
    SAFE_DELETE(this->frozenValues);
}


/*
 * view::View3D::Render
 */
void view::View3D::Render(void) {
    if (this->doHookCode()) {
        this->doBeforeRenderHook();
    }

    CallRender *cr = this->rendererSlot.CallAs<CallRender>();
    CallRender3D *cr3d = this->rendererSlot.CallAs<CallRender3D>();

    this->fpsCounter.FrameBegin();

    // clear viewport
    if (this->setViewport) {
        ::glViewport(0, 0,
            static_cast<GLsizei>(this->camParams->TileRect().Width()),
            static_cast<GLsizei>(this->camParams->TileRect().Height()));
    }

    if (this->backCol.IsDirty()) {
        this->backCol.ResetDirty();
        utility::ColourParser::FromString(this->backCol.Param<param::StringParam>()->Value(),
            this->bkgndCol[0], this->bkgndCol[1], this->bkgndCol[2]);
    }

    ::glClearColor(this->bkgndCol[0], this->bkgndCol[1],
        this->bkgndCol[2], 0.0f);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if ((cr == NULL) && (cr3d == NULL)) {
        return; // empty enought
    }

    // camera settings
    if (cr3d != NULL) {
        (*cr3d)(1); // GetExtents
        if (!(cr3d->AccessBoundingBoxes() == this->bboxs)) {
            this->bboxs = cr3d->AccessBoundingBoxes();

            if (this->firstImg) {
                this->ResetView();
                this->firstImg = false;
                if (!this->cameraSettingsSlot.Param<param::StringParam>()->Value().IsEmpty()) {
                    this->onRestoreCamera(this->restoreCameraSettingsSlot);
                }

            }
        }

        if (this->animPlay.Param<param::BoolParam>()->Value()) {
            float fps = this->animSpeed.Param<param::FloatParam>()->Value();
            float seconds = float(vislib::sys::PerformanceCounter::ToMillis(
                this->animTimer.Difference()) * 0.001);
            this->timeFrame += fps * seconds;
        }
        if (((unsigned int)this->timeFrame) >= cr3d->TimeFramesCount()) {
            this->timeFrame = 0.0f;
        }

        this->animTimer.SetMark();

        cr3d->SetTime(this->frozenValues ? this->frozenValues->time : this->timeFrame);
        cr3d->SetCameraParameters(this->cam.Parameters()); // < here we use the 'active' parameters!
        cr3d->SetLastFrameTime(this->fpsCounter.LastFrameTime());
    } else {
        cr->SetCameraParameters(this->cam.Parameters()); // < here we use the 'active' parameters!
    }
    this->camParams->CalcClipping(this->bboxs.ClipBox(), 0.1f);

    // set light parameters
    ::glEnable(GL_LIGHTING);    // TODO: check renderer capabilities
    ::glEnable(GL_LIGHT0);
    const float lp[4] = {
        -this->lightDir.X(),
        -this->lightDir.Y(),
        -this->lightDir.Z(),
        0.0f};
    const float la[4] = {0.1f, 0.1f, 0.1f, 1.0f}; // Could be configurable
    glLightfv(GL_LIGHT0, GL_AMBIENT, la);
    const float ld[4] = {0.9f, 0.9f, 0.9f, 1.0f}; // Could be configurable
    glLightfv(GL_LIGHT0, GL_DIFFUSE, ld);

    //printf("Cam%u\n", this->cam.Parameters()->SyncNumber());

    // setup matrices
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    this->cam.glMultProjectionMatrix();

    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();

    if (isCamLight) glLightfv(GL_LIGHT0, GL_POSITION, lp);

    this->cam.glMultViewMatrix();

    if (!isCamLight) glLightfv(GL_LIGHT0, GL_POSITION, lp);

    // render bounding box backside
    if (this->showBBox.Param<param::BoolParam>()->Value()) {
        this->renderBBoxBackside();
    }
    ::glPushMatrix();

    // call for render
    if (cr3d != NULL) {
        (*cr3d)(0);
    } else {
        (*cr)();
    }

    // render bounding box front
    ::glPopMatrix();
    if (this->showBBox.Param<param::BoolParam>()->Value()) {
        this->renderBBoxFrontside();
    }

    if (this->softCursor.Param<param::BoolParam>()->Value()) {
        this->renderSoftCursor();
    }

    if (true) {
        //::glMatrixMode(GL_PROJECTION);
        //::glLoadIdentity();
        //::glMatrixMode(GL_MODELVIEW);
        //::glLoadIdentity();
        unsigned int ticks = vislib::sys::GetTicksOfDay();
        if ((ticks < this->fpsOutputTimer) || (ticks >= this->fpsOutputTimer + 1000)) {
            this->fpsOutputTimer = ticks;
            printf("FPS: %f\n", this->fpsCounter.FPS());
        }
    }

    this->fpsCounter.FrameEnd();

}


/*
 * view::View3D::ResetView
 */
void view::View3D::ResetView(void) {
    using namespace vislib::graphics;
    VLTRACE(VISLIB_TRCELVL_INFO, "View3D::ResetView\n");

    this->camParams->SetClip(0.1f, 100.0f);
    this->camParams->SetApertureAngle(30.0f);
    this->camParams->SetProjection(
        vislib::graphics::CameraParameters::MONO_PERSPECTIVE);
    this->camParams->SetStereoParameters(0.3f, /* this is not so clear! */
        vislib::graphics::CameraParameters::LEFT_EYE,
        0.0f /* this is autofocus */);
    this->camParams->Limits()->LimitClippingDistances(0.01f, 0.1f);

    if (!this->bboxs.IsWorldSpaceBBoxValid()) {
        this->bboxs.SetWorldSpaceBBox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0);
    }
    float dist = (0.5f
        * sqrtf((this->bboxs.WorldSpaceBBox().Width() * this->bboxs.WorldSpaceBBox().Width())
        + (this->bboxs.WorldSpaceBBox().Depth() * this->bboxs.WorldSpaceBBox().Depth())
        + (this->bboxs.WorldSpaceBBox().Height() * this->bboxs.WorldSpaceBBox().Height())))
        / tanf(this->cam.Parameters()->HalfApertureAngle());
    SceneSpacePoint3D bbc = this->bboxs.WorldSpaceBBox().CalcCenter();

    this->camParams->SetView(
        bbc + SceneSpaceVector3D(0.0f, 0.0f, dist),
        bbc, SceneSpaceVector3D(0.0f, 1.0f, 0.0f));

    this->zoomer1.SetSpeed(dist);
}


/*
 * view::View3D::Resize
 */
void view::View3D::Resize(unsigned int width, unsigned int height) {
    this->camParams->SetVirtualViewSize(
        static_cast<vislib::graphics::ImageSpaceType>(width), 
        static_cast<vislib::graphics::ImageSpaceType>(height));
}


/*
 * view::View3D::SetInputModifier
 */
void view::View3D::SetInputModifier(mmcInputModifier mod, bool down) {
    unsigned int modId = 0;
    switch (mod) {
        case MMC_INMOD_SHIFT:
            modId = vislib::graphics::InputModifiers::MODIFIER_SHIFT;
            break;
        case MMC_INMOD_CTRL:
            modId = vislib::graphics::InputModifiers::MODIFIER_CTRL;
            break;
        case MMC_INMOD_ALT:
            modId = vislib::graphics::InputModifiers::MODIFIER_ALT;
            break;
        default: return;
    }
    this->modkeys.SetModifierState(modId, down);
}


/*
 * view::View3D::OnRenderView
 */
bool view::View3D::OnRenderView(Call& call) {
    bool redirtyColour = false;
    float oldBR = 0, oldBG = 0, oldBB = 0;
    view::CallRenderView *crv = dynamic_cast<view::CallRenderView *>(&call);
    if (crv == NULL) return false;

    this->setViewport = false;
    if (crv->IsProjectionSet() || crv->IsTileSet() || crv->IsViewportSet()) {
        this->cam.SetParameters(this->camOverrides);
        this->camOverrides.DynamicCast<ClusterCameraParamOverride>()
            ->SetOverrides(*crv);
    }
    if (crv->IsBackgroundSet()) {
        redirtyColour = this->backCol.IsDirty();
        this->backCol.ResetDirty();
        oldBR = this->bkgndCol[0];
        this->bkgndCol[0] = static_cast<float>(crv->BackgroundRed()) / 255.0f;
        oldBG = this->bkgndCol[1];
        this->bkgndCol[1] = static_cast<float>(crv->BackgroundGreen()) / 255.0f;
        oldBB = this->bkgndCol[2];
        this->bkgndCol[2] = static_cast<float>(crv->BackgroundBlue()) / 255.0f;
    }

    this->Render();

    if (crv->IsProjectionSet() || crv->IsTileSet() || crv->IsViewportSet()) {
        this->cam.SetParameters(this->frozenValues ? this->frozenValues->camParams : this->camParams);
    }
    if (crv->IsBackgroundSet()) {
        if (redirtyColour) {
            this->backCol.ForceSetDirty();
        }
        this->bkgndCol[0] = oldBR;
        this->bkgndCol[1] = oldBG;
        this->bkgndCol[2] = oldBB;
    }
    this->setViewport = true;

    return true;
}


/*
 * view::View3D::UpdateFreeze
 */
void view::View3D::UpdateFreeze(bool freeze) {
    //printf("%s view\n", freeze ? "Freezing" : "Unfreezing");
    if (freeze) {
        if (this->frozenValues == NULL) {
            this->frozenValues = new FrozenValues();
            this->camOverrides.DynamicCast<ClusterCameraParamOverride>()
                ->SetParametersBase(this->frozenValues->camParams);
            this->cam.SetParameters(this->frozenValues->camParams);
        }
        *(this->frozenValues->camParams) = *this->camParams;
        this->frozenValues->time = this->timeFrame;
    } else {
        this->camOverrides.DynamicCast<ClusterCameraParamOverride>()
            ->SetParametersBase(this->camParams);
        this->cam.SetParameters(this->camParams);
        SAFE_DELETE(this->frozenValues);
    }
}


/*
 * view::View3D::create
 */
bool view::View3D::create(void) {
    
    this->cursor2d.SetButtonCount(3); /* This could be configurable. */
    this->modkeys.SetModifierCount(3);

    this->rotator1.SetCameraParams(this->camParams);
    this->rotator1.SetTestButton(0 /* left mouse button */);
    this->rotator1.SetAltModifier(
        vislib::graphics::InputModifiers::MODIFIER_SHIFT);
    this->rotator1.SetModifierTestCount(1);
    this->rotator1.SetModifierTest(0,
        vislib::graphics::InputModifiers::MODIFIER_CTRL, false);

    this->rotator2.SetCameraParams(this->camParams);
    this->rotator2.SetTestButton(0 /* left mouse button */);
    this->rotator2.SetAltModifier(
        vislib::graphics::InputModifiers::MODIFIER_SHIFT);
    this->rotator2.SetModifierTestCount(1);
    this->rotator2.SetModifierTest(0,
        vislib::graphics::InputModifiers::MODIFIER_CTRL, true);

    this->zoomer1.SetCameraParams(this->camParams);
    this->zoomer1.SetTestButton(2 /* mid mouse button */);
    this->zoomer1.SetModifierTestCount(1);
    this->zoomer1.SetModifierTest(0,
        vislib::graphics::InputModifiers::MODIFIER_CTRL, false);

    this->zoomer2.SetCameraParams(this->camParams);
    this->zoomer2.SetTestButton(2 /* mid mouse button */);
    this->zoomer2.SetModifierTestCount(1);
    this->zoomer2.SetModifierTest(0,
        vislib::graphics::InputModifiers::MODIFIER_CTRL, true);

    this->cursor2d.SetCameraParams(this->camParams);
    this->cursor2d.RegisterCursorEvent(&this->rotator1);
    this->cursor2d.RegisterCursorEvent(&this->rotator2);
    this->cursor2d.RegisterCursorEvent(&this->zoomer1);
    this->cursor2d.RegisterCursorEvent(&this->zoomer2);
    this->cursor2d.SetInputModifiers(&this->modkeys);

    this->modkeys.RegisterObserver(&this->cursor2d);

    this->firstImg = true;

    return true;
}


/*
 * view::View3D::release
 */
void view::View3D::release(void) {
    this->cursor2d.UnregisterCursorEvent(&this->rotator1);
    this->cursor2d.UnregisterCursorEvent(&this->rotator2);
    this->cursor2d.UnregisterCursorEvent(&this->zoomer1);
    this->cursor2d.UnregisterCursorEvent(&this->zoomer2);
    SAFE_DELETE(this->frozenValues);
}


/*
 * view::View3D::renderBBox
 */
void view::View3D::renderBBox(void) {
    const vislib::math::Cuboid<float>& boundingBox
        = this->bboxs.WorldSpaceBBox();
    if (!this->bboxs.IsWorldSpaceBBoxValid()) {
        this->bboxs.SetWorldSpaceBBox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f);
    }

    ::glBegin(GL_QUADS);

    ::glEdgeFlag(true);

    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Back());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Back());

    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Front());
    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Front());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Front());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Front());

    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Back());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Front());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Front());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Back());

    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Front());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Front());

    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Back());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Bottom(), boundingBox.Front());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Front());
    ::glVertex3f(boundingBox.Left(),  boundingBox.Top(),    boundingBox.Back());

    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Back());
    ::glVertex3f(boundingBox.Right(), boundingBox.Top(),    boundingBox.Front());
    ::glVertex3f(boundingBox.Right(), boundingBox.Bottom(), boundingBox.Front());

    ::glEnd();
}


/*
 * view::View3D::renderBBoxBackside
 */
void view::View3D::renderBBoxBackside(void) {
    ::glDisable(GL_LIGHTING);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_CULL_FACE);
    ::glCullFace(GL_FRONT);
    ::glEnable(GL_LINE_SMOOTH);
    ::glLineWidth(1.25f);
    ::glDisable(GL_TEXTURE_2D);
    ::glPolygonMode(GL_BACK, GL_LINE);

    ::glColor4ub(255, 255, 255, 50);
    this->renderBBox();

    ::glPolygonMode(GL_BACK, GL_FILL);
    ::glDisable(GL_DEPTH_TEST);

    ::glColor4ub(255, 255, 255, 12);
    this->renderBBox();

    ::glCullFace(GL_BACK);
}


/*
 * view::View3D::renderBBoxFrontside
 */
void view::View3D::renderBBoxFrontside(void) {
    ::glDisable(GL_LIGHTING);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    ::glEnable(GL_DEPTH_TEST);
    ::glDepthFunc(GL_LEQUAL);
    ::glEnable(GL_CULL_FACE);
    ::glCullFace(GL_BACK);
    ::glEnable(GL_LINE_SMOOTH);
    ::glLineWidth(1.75f);
    ::glDisable(GL_TEXTURE_2D);
    ::glPolygonMode(GL_FRONT, GL_LINE);

    ::glColor4ub(255, 255, 255, 100);
    this->renderBBox();

    ::glDepthFunc(GL_LESS);
    ::glPolygonMode(GL_FRONT, GL_FILL);
}


/*
 * view::View3D::onSetCursor2DButtonState
 */
bool view::View3D::onSetCursor2DButtonState(Call& call) {
    view::CallCursorInput *cci = dynamic_cast<view::CallCursorInput *>(&call);
    if (cci == NULL) return false;
    this->SetCursor2DButtonState(cci->Btn(), cci->Down());
    return true;
}


/*
 * view::View3D::onSetCursor2DPosition
 */
bool view::View3D::onSetCursor2DPosition(Call& call) {
    view::CallCursorInput *cci = dynamic_cast<view::CallCursorInput *>(&call);
    if (cci == NULL) return false;

    float w = this->camParams->VirtualViewSize().Width();
    float h = this->camParams->VirtualViewSize().Height();
    float s = vislib::math::Min(w, h);
    this->SetCursor2DPosition(0.5f * w + cci->X() * s, 0.5f * h + cci->Y() * s);

    return true;
}


/*
 * view::View3D::onSetInputModifier
 */
bool view::View3D::onSetInputModifier(Call& call) {
    view::CallCursorInput *cci = dynamic_cast<view::CallCursorInput *>(&call);
    if (cci == NULL) return false;
    this->SetInputModifier(cci->Mod(), cci->Down());
    return true;
}


/*
 * view::View3D::onResetView
 */
bool view::View3D::onResetView(Call& call) {
    this->ResetView();
    return true;
}


/*
 * view::View3D::renderSoftCursor
 */
void view::View3D::renderSoftCursor(void) {
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();

    // transform to tile coordinates!
    ::glTranslatef(-1.0, -1.0f, -1.0f);
    ::glScalef(2.0f, 2.0f, 1.0f);
    ::glScalef(1.0f / this->cam.Parameters()->TileRect().Width(),
        1.0f / this->cam.Parameters()->TileRect().Height(), 1.0f);
    ::glTranslatef(-this->cam.Parameters()->TileRect().Left(),
        -this->cam.Parameters()->TileRect().Bottom(), 0.0f);

    float w = this->camParams->VirtualViewSize().Width();
    float h = this->camParams->VirtualViewSize().Height();
    float s = vislib::math::Min(w, h);
    float x = this->cursor2d.X();
    float y = this->cursor2d.Y();

    x = (x - 0.5f * w) / s;
    y = (y - 0.5f * h) / s;

    w = this->cam.Parameters()->VirtualViewSize().Width();
    h = this->cam.Parameters()->VirtualViewSize().Height();
    s = vislib::math::Min(w, h);

    x = x * s + 0.5f * w;
    y = y * s + 0.5f * h;

    ::glTranslatef(x, y, 0.0f);
    ::glScalef(1.5f, -1.5f, 1.0f);

    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glEnable(GL_LINE_SMOOTH);
    ::glLineWidth(1.5f);
    ::glBindTexture(GL_TEXTURE_2D, 0);
    ::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glDisable(GL_DEPTH_TEST);

    ::glBegin(GL_TRIANGLE_FAN);
    ::glColor3ub(255, 255, 255); ::glVertex2i(0, 0);
    ::glColor3ub(245, 245, 245); ::glVertex2i(0, 17);
    ::glColor3ub(238, 238, 238); ::glVertex2i(4, 13);
    ::glColor3ub(211, 211, 211); ::glVertex2i(7, 18);
    ::glVertex2i(9, 18);
    ::glVertex2i(9, 16);
    ::glColor3ub(234, 234, 234); ::glVertex2i(7, 12);
    ::glColor3ub(226, 226, 226); ::glVertex2i(12, 12);
    ::glEnd();
    ::glBegin(GL_LINE_LOOP);
    ::glColor3ub(0, 0, 0);
    ::glVertex2i(0, 0);
    ::glVertex2i(0, 17);
    ::glVertex2i(4, 13);
    ::glVertex2i(7, 18);
    ::glVertex2i(9, 18);
    ::glVertex2i(9, 16);
    ::glVertex2i(7, 12);
    ::glVertex2i(12, 12);
    ::glEnd();

    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable(GL_BLEND);
}


/*
 * view::View3D::onStoreCamera
 */
bool view::View3D::onStoreCamera(param::ParamSlot& p) {
    vislib::TStringSerialiser strser;
    strser.ClearData();
    this->camParams->Serialise(strser);
    vislib::TString str(strser.GetString());
    str.EscapeCharacters(_T('\\'), _T("\n\r\t"), _T("nrt"));
    str.Append(_T("}"));
    str.Prepend(_T("{"));
    this->cameraSettingsSlot.Param<param::StringParam>()->SetValue(str);

    vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_INFO,
        "Camera parameters stored in \"%s\"",
        this->cameraSettingsSlot.FullName().PeekBuffer());
    return true;
}


/*
 * view::View3D::onRestoreCamera
 */
bool view::View3D::onRestoreCamera(param::ParamSlot& p) {
    vislib::TString str = this->cameraSettingsSlot.Param<param::StringParam>()->Value();
    try {
        if ((str[0] != _T('{')) || (str[str.Length() - 1] != _T('}'))) {
            throw vislib::Exception("invalid string: surrounding brackets missing", __FILE__, __LINE__);
        }
        str = str.Substring(1, str.Length() - 2);
        if (!str.UnescapeCharacters(_T('\\'), _T("\n\r\t"), _T("nrt"))) {
            throw vislib::Exception("unrecognised escape sequence", __FILE__, __LINE__);
        }
        vislib::TStringSerialiser strser(str);
        vislib::graphics::CameraParamsStore cps;
        cps = *this->camParams.operator->();
        cps.Deserialise(strser);
        cps.SetVirtualViewSize(this->camParams->VirtualViewSize());
        cps.SetTileRect(this->camParams->TileRect());
        *this->camParams.operator->() = cps;
        // now avoid resetting the camera by the initialisation
        if (!this->bboxs.IsWorldSpaceBBoxValid()) {
            this->bboxs.SetWorldSpaceBBox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f);
        }

        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_INFO,
            "Camera parameters restored from \"%s\"",
            this->cameraSettingsSlot.FullName().PeekBuffer());
    } catch(vislib::Exception ex) {
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
            "Cannot restore camera parameters from \"%s\": %s (%s; %d)",
            this->cameraSettingsSlot.FullName().PeekBuffer(),
            ex.GetMsgA(), ex.GetFile(), ex.GetLine());
    } catch(...) {
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
            "Cannot restore camera parameters from \"%s\": unexpected exception",
            this->cameraSettingsSlot.FullName().PeekBuffer());
    }
    return true;
}


/*
 * view::View3D::onResetView
 */
bool view::View3D::onResetView(param::ParamSlot& p) {
    this->ResetView();
    return true;
}
