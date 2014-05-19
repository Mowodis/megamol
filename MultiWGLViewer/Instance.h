/*
 * Instance.h
 *
 * Copyright (C) 2006 - 2011 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_WGL_INSTANCE_H_INCLUDED
#define MEGAMOL_WGL_INSTANCE_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "ApiHandle.h"


namespace megamol {
namespace wgl {

    /** 
     * Library instance
     */
    class Instance : public ApiHandle {
    public:

        /** The window class name */
        static const TCHAR* WindowClassName;

        /**
         * Gets the application instance handle
         *
         * @return The application instance handle
         */
        static HINSTANCE HInstance(void) {
            return hInst;
        }

        /** Ctor */
        Instance(void);

        /** Dtor */
        virtual ~Instance(void);

        /**
         * Initializes the instance
         *
         * @param hInst The application instance handle
         *
         * @return True on success
         */
        bool Init(HINSTANCE hInst);

        /**
         * Processes Window event
         *
         * @return True if the application should continue
         */
        bool ProcessEvents(void);

    private:

        /** The application instance handle */
        static HINSTANCE hInst;

        /** Reference counter how many instances use the window class */
        static unsigned int wndClsRegistered;

        /**
         * Deactivates Desktop Window Composition
         */
        void deactivateCompositeDesktop(void);

        /** Flag whether or not the application is running */
        bool running;

    };


} /* end namespace wgl */
} /* end namespace megamol */

#endif /* MEGAMOL_WGL_INSTANCE_H_INCLUDED */
