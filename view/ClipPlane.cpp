/*
 * ClipPlane.cpp
 *
 * Copyright (C) 2009 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "ClipPlane.h"
#include "param/BoolParam.h"
#include "param/StringParam.h"
#include "param/Vector3fParam.h"
#include "utility/ColourParser.h"
#include "view/CallClipPlane.h"
#include "vislib/Log.h"
#include "vislib/mathfunctions.h"
#include "vislib/ShallowPoint.h"
#include "vislib/Vector.h"

using namespace megamol::core;


/*
 * view::ClipPlane::ClipPlane
 */
view::ClipPlane::ClipPlane(void) : Module(),
        getClipPlaneSlot("getclipplane", "Provides the clipping plane"),
        plane(),
        enableSlot("enable", "Disables or enables the clipping plane"),
        colourSlot("colour", "Defines the colour of the clipping plane"),
        normalSlot("normal", "Defines the normal of the clipping plane"),
        pointSlot("point", "Defines a point in the clipping plane") {

    this->plane.Set(vislib::math::Point<float, 3>(0.0, 0.0f, 0.0f),
        vislib::math::Vector<float, 3>(1.0f, 0.0f, 0.0f));
    this->col[0] = this->col[1] = this->col[2] = 128;

    view::CallClipPlaneDescription ccpd;
    this->getClipPlaneSlot.SetCallback(ccpd.ClassName(), ccpd.FunctionName(0),
        &ClipPlane::requestPlane);
    this->MakeSlotAvailable(&this->getClipPlaneSlot);

    this->enableSlot << new param::BoolParam(false);
    this->MakeSlotAvailable(&this->enableSlot);

    this->colourSlot << new param::StringParam(
        utility::ColourParser::ToString(
            static_cast<float>(this->col[0]) / 255.0f,
            static_cast<float>(this->col[1]) / 255.0f,
            static_cast<float>(this->col[2]) / 255.0f));
    this->MakeSlotAvailable(&this->colourSlot);

    this->normalSlot << new param::Vector3fParam(this->plane.Normal());
    this->MakeSlotAvailable(&this->normalSlot);

    this->pointSlot << new param::Vector3fParam(
        vislib::math::Vector<float, 3>(this->plane.Point()));
    this->MakeSlotAvailable(&this->pointSlot);
}


/*
 * view::ClipPlane::~ClipPlane
 */
view::ClipPlane::~ClipPlane(void) {
    this->Release();
}


/*
 * view::ClipPlane::create
 */
bool view::ClipPlane::create(void) {
    // intentionally empty
    return true;
}


/*
 * view::ClipPlane::release
 */
void view::ClipPlane::release(void) {
    // intentionally empty
}


/*
 * view::ClipPlane::requestPlane
 */
bool view::ClipPlane::requestPlane(Call& call) {
    view::CallClipPlane *ccp = dynamic_cast<view::CallClipPlane*>(&call);
    if (ccp == NULL) return false;

    if (!this->enableSlot.Param<param::BoolParam>()->Value()) {
        // clipping plane is disabled
        return false;
    }

    if (this->colourSlot.IsDirty()) {
        this->colourSlot.ResetDirty();
        float r, g, b;
        if (utility::ColourParser::FromString(this->colourSlot.Param<param::StringParam>()->Value(), r, g, b)) {
            this->col[0] = static_cast<unsigned char>(vislib::math::Clamp(r, 0.0f, 1.0f) * 255.0f);
            this->col[1] = static_cast<unsigned char>(vislib::math::Clamp(g, 0.0f, 1.0f) * 255.0f);
            this->col[2] = static_cast<unsigned char>(vislib::math::Clamp(b, 0.0f, 1.0f) * 255.0f);
        } else {
            vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
                "Unable to parse colour");
        }
    }

    if (this->normalSlot.IsDirty() || this->pointSlot.IsDirty()) {
        this->normalSlot.ResetDirty();
        this->pointSlot.ResetDirty();
        this->plane.Set(vislib::math::ShallowPoint<float, 3>(
                const_cast<float*>(
                this->pointSlot.Param<param::Vector3fParam>()->Value().PeekComponents())),
            this->normalSlot.Param<param::Vector3fParam>()->Value());
    }

    ccp->SetColour(this->col[0], this->col[1], this->col[2]);
    ccp->SetPlane(this->plane);

    return true;
}
