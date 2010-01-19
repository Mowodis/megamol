/*
 * ModuleDescriptionManager.h
 *
 * Copyright (C) 2008 by Universitaet Stuttgart (VIS). 
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_MODULEDESCRIPTIONMANAGER_H_INCLUDED
#define MEGAMOLCORE_MODULEDESCRIPTIONMANAGER_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */


#include "ObjectDescriptionManager.h"
#include "ModuleDescription.h"


namespace megamol {
namespace core {


    /**
     * Class of rendering graph module description manager
     */
    class ModuleDescriptionManager
        : public ObjectDescriptionManager<ModuleDescription> {
    public:

        /**
         * Returns the only instance of this class.
         *
         * @return The only instance of this class.
         */
        static ModuleDescriptionManager * Instance();

    private:

        /** Private ctor. */
        ModuleDescriptionManager(void);

        /** Private dtor. */
        virtual ~ModuleDescriptionManager(void);

        /**
         * Registers a module (if available) description. Template parameter
         * Cp is the module description class.
         */
        template<class Cp>
        void registerDescription(void);

        /**
         * Registers a module (if available) description.
         *
         * @param desc The module description to be registered. The memory of
         *             the description object will be managed by the manager
         *             object. The caller must not alter the memory (e.g. free)
         *             after this call was issued.
         */
        void registerDescription(ModuleDescription* desc);

        /**
         * Registers a module (if available) using an auto-description object.
         * Template parameter Cp is the module class to be described by an
         * auto-description object.
         */
        template<class Cp>
        void registerAutoDescription(void);

    };


} /* end namespace core */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_MODULEDESCRIPTIONMANAGER_H_INCLUDED */
