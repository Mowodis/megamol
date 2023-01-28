/*
 * PeerDisconnectedException.h
 *
 * Copyright (C) 2006 - 2010 by Visualisierungsinstitut Universitaet Stuttgart.
 * Alle Rechte vorbehalten.
 */

#ifndef VISLIB_PEERDISCONNECTEDEXCEPTION_H_INCLUDED
#define VISLIB_PEERDISCONNECTEDEXCEPTION_H_INCLUDED
#pragma once
#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(push, off)
#endif /* defined(_WIN32) && defined(_MANAGED) */


#include "vislib/Exception.h"
#include "vislib/String.h"


namespace vislib::net {


/**
 * This exception is thrown if an AbstractCommChannel or a derived class
 * detects a graceful disconnect of the peer node.
 */
class PeerDisconnectedException : public Exception {

public:
    /**
     * Creates an exception message for the local end point description
     * 'localEndPoint'.
     *
     * @param localEndPoint A description for the local end point which
     *                      was disconnected. It is recommended using the
     *                      bind address for this.
     *
     * @return The default message to be used for a
     *         PeerDisconnectedException.
     */
    static StringA FormatMessageForLocalEndpoint(const char* localEndPoint);

    /**
     * Creates an exception message for the local end point description
     * 'localEndPoint'.
     *
     * @param localEndPoint A description for the local end point which
     *                      was disconnected. It is recommended using the
     *                      bind address for this.
     *
     * @return The default message to be used for a
     *         PeerDisconnectedException.
     */
    static StringW FormatMessageForLocalEndpoint(const wchar_t* localEndPoint);

    /**
     * Create a new exception. The ownership of the memory designated by
     * 'msg' and 'file' remains at the caller, the class creates a deep
     * copy.
     *
     * @param msg  A description of the exception.
     * @param file The file the exception was thrown in.
     * @param line The line the exception was thrown in.
     */
    PeerDisconnectedException(const char* msg, const char* file, const int line);

    /**
     * Create a new exception. The ownership of the memory designated by
     * 'msg' and 'file' remains at the caller, the class creates a deep
     * copy.
     *
     * @param msg  A description of the exception.
     * @param file The file the exception was thrown in.
     * @param line The line the exception was thrown in.
     */
    PeerDisconnectedException(const wchar_t* msg, const char* file, const int line);

    /**
     * Create a clone of 'rhs'.
     *
     * @param rhs The object to be cloned.
     */
    PeerDisconnectedException(const PeerDisconnectedException& rhs);

    /** Dtor. */
    ~PeerDisconnectedException() override;

private:
    /** Superclass typedef. */
    typedef Exception Super;
};

} // namespace vislib::net

#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(pop)
#endif /* defined(_WIN32) && defined(_MANAGED) */
#endif /* VISLIB_PEERDISCONNECTEDEXCEPTION_H_INCLUDED */
