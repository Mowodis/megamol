/*
 * ShaderSourceFactory.cpp
 *
 * Copyright (C) 2008 by Universitaet Stuttgart (VISUS). 
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "ShaderSourceFactory.h"
#include "utility/BTFParser.h"
#include "utility/xml/XmlReader.h"
#include "vislib/StringTokeniser.h"
#include "vislib/Log.h"
#include "vislib/Path.h"

using namespace megamol::core;


/*
 * utility::ShaderSourceFactory::ShaderSourceFactory
 */
utility::ShaderSourceFactory::ShaderSourceFactory(utility::Configuration& config)
        : config(config), root(), fileIds() {
    // intentionally empty
}


/*
 * utility::ShaderSourceFactory::~ShaderSourceFactory
 */
utility::ShaderSourceFactory::~ShaderSourceFactory(void) {
    // intentionally empty (ATM)
}


/*
 * utility::ShaderSourceFactory::LoadBTF
 */
bool utility::ShaderSourceFactory::LoadBTF(const vislib::StringA & name) {
    using vislib::sys::Log;
    if (name.Find("::") != vislib::StringA::INVALID_POS) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "BTF root namespace \"%s\" is illegal: contains \"::\"",
            name.PeekBuffer());
        return false;
    }
    if (!this->root.FindChild(name).IsNull()) {
        //Log::DefaultLog.WriteMsg(Log::LEVEL_INFO + 750,
        //    "BTF \"%s\" already loaded.", name.PeekBuffer());
        return true;
    }

    vislib::StringW path;
    mmcValueType type;
    const void *data = this->config.GetValue(MMC_CFGID_SHADER_DIR,
        (const wchar_t*)NULL, &type);
    if (type == MMC_TYPE_CSTR) {
        path = vislib::StringA(static_cast<const char *>(data));
    } else if (type == MMC_TYPE_WSTR) {
        path = vislib::StringW(static_cast<const wchar_t *>(data));
    } else {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load btf \"%s\": shader directory unavailable\n",
            name.PeekBuffer());
        return false;
    }
    vislib::StringW filename = vislib::sys::Path::Concatenate(path,
        vislib::StringW(name));
    filename.Append(L".btf");

    if (!vislib::sys::File::Exists(filename)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load btf \"%s\": file not found\n",
            name.PeekBuffer());
        return false;
    }

#if defined(_WIN32) && (defined(DEBUG) || defined(_DEBUG))
    WIN32_FIND_DATAW findData;
    HANDLE hFind = ::FindFirstFileW(filename.PeekBuffer(), &findData);
    if (hFind != INVALID_HANDLE_VALUE) {
        vislib::StringW osfilename(findData.cFileName);
        ::FindClose(hFind);
        osfilename.Truncate(osfilename.Length() - 4);
        if (!osfilename.Equals(name)) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
                "Unable to load btf \"%s\": wrong file name case\n",
                name.PeekBuffer());
            return false;
        }
    } else {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load btf \"%s\": file not found(2)\n",
            name.PeekBuffer());
        return false;
    }
#endif /* _WIN32 && (DEBUG || _DEBUG) */

    xml::XmlReader reader;
    if (!reader.OpenFile(filename)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load btf \"%s\": cannot open file\n",
            name.PeekBuffer());
        return false;
    }

    vislib::SmartPtr<BTFParser::BTFElement> fileroot = new BTFParser::BTFNamespace(name);
    this->root.Children().Add(fileroot);

    BTFParser parser(*this);
    parser.SetRootElement(fileroot, new ShallowBTFNamespace(this->root));

    if (!parser.Parse(reader)) {

        vislib::SingleLinkedList<vislib::StringA>::Iterator mi = parser.Messages();
        while (mi.HasNext()) {
            Log::DefaultLog.WriteMsg(Log::LEVEL_INFO, "Parser: %s",
                mi.Next().PeekBuffer());
        }

        this->root.Children().RemoveAll(fileroot);
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load btf \"%s\": failed to parse\n",
            name.PeekBuffer());
        return false;
    }

    vislib::SingleLinkedList<vislib::StringA>::Iterator mi = parser.Messages();
    while (mi.HasNext()) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_INFO, "Parser: %s",
            mi.Next().PeekBuffer());
    }
    reader.CloseFile();

    if (this->fileIds.Find(name) == NULL) {
        this->fileIds.Append(name);
        Log::DefaultLog.WriteMsg(Log::LEVEL_INFO,
            "BTF \"%s\" indexed as %d\n",
            name.PeekBuffer(), this->fileIds.Count());
    }

    return true;
}


/*
 * utility::ShaderSourceFactory::MakeShaderSnippet
 */
vislib::SmartPtr<vislib::graphics::gl::ShaderSource::Snippet>
utility::ShaderSourceFactory::MakeShaderSnippet(const vislib::StringA& name) {
    vislib::SmartPtr<BTFParser::BTFElement> el = this->getBTFElement(name);
    if (el.IsNull()) return NULL;
    return this->makeSnippet(el.DynamicCast<BTFParser::BTFSnippet>());
}


/*
 * utility::ShaderSourceFactory::MakeShaderSource
 */
bool utility::ShaderSourceFactory::MakeShaderSource(const vislib::StringA& name,
        vislib::graphics::gl::ShaderSource& outShaderSrc) {
    using vislib::sys::Log;
    vislib::SmartPtr<BTFParser::BTFElement> el = this->getBTFElement(name);
    BTFParser::BTFShader *s = el.DynamicCast<BTFParser::BTFShader>();
    outShaderSrc.Clear();
    if (s == NULL) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to find requested shader \"%s\"\n", name.PeekBuffer());
        return false;
    }
    return this->makeShaderSource(s, outShaderSrc);
}


/*
 * utility::ShaderSourceFactory::fileNo
 */
int utility::ShaderSourceFactory::fileNo(const vislib::StringA& name) {
    vislib::StringA::Size pos = name.Find("::");
    INT_PTR idx = 0;
    if (pos == vislib::StringA::INVALID_POS) {
        idx = this->fileIds.IndexOf(name);
    } else {
        idx = this->fileIds.IndexOf(name.Substring(0, pos));
    }
    if (idx == vislib::Array<vislib::StringA>::INVALID_POS) {
        idx = 0;
    } else {
        idx++;
    }
    return static_cast<int>(idx);
}


/*
 * utility::ShaderSourceFactory::getBTFElement
 */
vislib::SmartPtr<utility::BTFParser::BTFElement>
utility::ShaderSourceFactory::getBTFElement(const vislib::StringA& name) {
    using vislib::sys::Log;
    vislib::Array<vislib::StringA> namepath = vislib::StringTokeniserA::Split(name, "::");
    if ((namepath.Count() > 0) && namepath[0].IsEmpty()) {
        namepath.RemoveFirst(); // allow operator global namespace
    }
    if (namepath.Count() < 2) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "Unable to load shader snippet \"%s\": incomplete name path\n",
            name.PeekBuffer());
        return NULL;
    }
    if (!this->LoadBTF(namepath[0])) return NULL;
    vislib::SmartPtr<BTFParser::BTFElement> e = this->root.FindChild(namepath[0]);
    for (SIZE_T i = 1; i < namepath.Count(); i++) {
        if (e.IsNull()) return NULL;
        if (e.DynamicCast<BTFParser::BTFNamespace>() == NULL) return NULL;
        e = e.DynamicCast<BTFParser::BTFNamespace>()->FindChild(namepath[i]);
    }
    return e;
}


/*
 * utility::ShaderSourceFactory::makeShaderSource
 */
bool utility::ShaderSourceFactory::makeShaderSource(
        utility::BTFParser::BTFShader *s,
        vislib::graphics::gl::ShaderSource& o) {
    ASSERT(s != NULL);

    unsigned int idx = 0;
    vislib::SingleLinkedList<vislib::SmartPtr<BTFParser::BTFElement>
        >::Iterator i = s->Children().GetIterator();
    while (i.HasNext()) {
        vislib::SmartPtr<BTFParser::BTFElement> e = i.Next();
        if (e.IsNull()) {
            vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_WARN,
                "Ignoring empty shader element.");
            idx++;
            continue;
        }

        BTFParser::BTFSnippet *c = e.DynamicCast<BTFParser::BTFSnippet>();
        if (c != NULL) {
            vislib::SingleLinkedList<vislib::StringA> keys = s->NameIDs().FindKeys(idx);
            vislib::SmartPtr<vislib::graphics::gl::ShaderSource::Snippet>
                snip = this->makeSnippet(c);
            if (snip.IsNull()) {
                return false;
            }
            if (keys.IsEmpty()) {
                o.Append(snip);
            } else {
                if (keys.Count() > 1) {
                    vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_WARN,
                        "Multiple shader snippet names encountered. Choosing \"%s\"",
                        keys.First().PeekBuffer());
                }
                o.Append(keys.First(), snip);
            }
        } else {
            BTFParser::BTFShader *sh = e.DynamicCast<BTFParser::BTFShader>();
            if (sh != NULL) {
                if (!this->makeShaderSource(sh, o)) {
                    return false;
                }
            } else {
                vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_WARN,
                    "Ignoring unknown shader element \"%s\".",
                    e->Name().PeekBuffer());
            }
        }

        idx++;
    }
    return true;
}


/*
 * utility::ShaderSourceFactory::makeSnippet
 */
vislib::SmartPtr<vislib::graphics::gl::ShaderSource::Snippet>
utility::ShaderSourceFactory::makeSnippet(utility::BTFParser::BTFSnippet *s) {
    if (s == NULL) return NULL;

    if (s->Type() == BTFParser::BTFSnippet::VERSION_SNIPPET) {
        try {
            return new vislib::graphics::gl::ShaderSource::VersionSnippet(
                vislib::CharTraitsA::ParseInt(s->Content()));
        } catch(...) {
            vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_ERROR,
                "Unable to create version shader code snippet.\n");
        }
    }

    if (s->Type() == BTFParser::BTFSnippet::EMPTY_SNIPPET) {
        vislib::sys::Log::DefaultLog.WriteMsg(vislib::sys::Log::LEVEL_WARN,
            "Empty shader code snippet encountered.\n");
    }

    vislib::StringA content;
    content.Format("\n#line %d %d\n%s", int(s->Line() - 1),
        this->fileNo(s->File()), s->Content().PeekBuffer());

    return new vislib::graphics::gl::ShaderSource::StringSnippet(content);
}
