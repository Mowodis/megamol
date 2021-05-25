/*
 * GUIUtils.h
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */


#ifndef MEGAMOL_GUI_GUIUTILS_INCLUDED
#define MEGAMOL_GUI_GUIUTILS_INCLUDED
#pragma once


#define IMGUI_DEFINE_MATH_OPERATORS
#define IMGUI_DISABLE_OBSOLETE_FUNCTIONS
#include "imgui.h"
#include "imgui_internal.h"

#include "mmcore/param/AbstractParamPresentation.h"
#include "mmcore/utility/log/Log.h"
#include "mmcore/view/Input.h"
#include "vislib/UTF8Encoder.h"
#include "vislib/math/Ternary.h"


/// #define GUI_VERBOSE


#define GUI_INVALID_ID (UINT_MAX)
#define GUI_SLOT_RADIUS (8.0f * megamol::gui::gui_scaling.Get())
#define GUI_LINE_THICKNESS (3.0f * megamol::gui::gui_scaling.Get())
#define GUI_RECT_CORNER_RADIUS (0.0f)
#define GUI_MAX_MULITLINE (8)
#define GUI_DND_CALLSLOT_UID_TYPE ("DND_CALL")
#define GUI_GRAPH_BORDER (GUI_SLOT_RADIUS * 4.0f)
#define GUI_MULTISELECT_MODIFIER (ImGui::GetIO().KeyShift)

#define GUI_JSON_TAG_GUI ("GUIState")
#define GUI_JSON_TAG_WINDOW_CONFIGS ("WindowConfigurations")
#define GUI_JSON_TAG_CONFIGURATOR ("ConfiguratorState")
#define GUI_JSON_TAG_GRAPHS ("GraphStates")
#define GUI_JSON_TAG_PROJECT ("Project")
#define GUI_JSON_TAG_MODULES ("Modules")
#define GUI_JSON_TAG_INTERFACES ("Interfaces")
/// GUI_JSON_TAG_GUISTATE_PARAMETERS ("ParameterStates") is defined in AbstractParamPresentation.h

#define GUI_START_TAG_SET_GUI_STATE ("mmSetGUIState([=[")
#define GUI_END_TAG_SET_GUI_STATE ("]=])")
#define GUI_START_TAG_SET_GUI_VISIBILITY ("mmSetGUIVisible(")
#define GUI_END_TAG_SET_GUI_VISIBILITY (")")
#define GUI_START_TAG_SET_GUI_SCALE ("mmSetGUIScale(")
#define GUI_END_TAG_SET_GUI_SCALE (")")

// Global Colors
#define GUI_COLOR_TEXT_ERROR (ImVec4(0.9f, 0.1f, 0.0f, 1.0f))
#define GUI_COLOR_TEXT_WARN (ImVec4(0.75f, 0.75f, 0.0f, 1.0f))

#define GUI_COLOR_BUTTON_MODIFIED (ImVec4(0.75f, 0.0f, 0.25f, 1.0f))
#define GUI_COLOR_BUTTON_MODIFIED_HIGHLIGHT (ImVec4(0.9f, 0.0f, 0.25f, 1.0f))

#define GUI_COLOR_SLOT_CALLER (ImVec4(0.0f, 0.75f, 1.0f, 1.0f))
#define GUI_COLOR_SLOT_CALLEE (ImVec4(0.75f, 0.0f, 1.0f, 1.0f))
#define GUI_COLOR_SLOT_COMPATIBLE (ImVec4(0.5f, 0.9f, 0.25f, 1.0f))
#define GUI_COLOR_SLOT_REQUIRED (ImVec4(0.9f, 0.75f, 0.1f, 1.0f))

#define GUI_COLOR_GROUP_HEADER (ImVec4(0.0f, 0.5f, 0.25f, 1.0f))
#define GUI_COLOR_GROUP_HEADER_HIGHLIGHT (ImVec4(0.0f, 0.75f, 0.5f, 1.0f))

// Texture File Names
#define GUI_DEFAULT_FONT_ROBOTOSANS ("Roboto-Regular.ttf")
#define GUI_DEFAULT_FONT_SOURCECODEPRO ("SourceCodePro-Regular.ttf")
#define GUI_TRANSPORT_ICON_PLAY ("transport_ctrl_play.png")
#define GUI_TRANSPORT_ICON_PAUSE ("transport_ctrl_pause.png")
#define GUI_TRANSPORT_ICON_FAST_FORWARD ("transport_ctrl_fast-forward.png")
#define GUI_TRANSPORT_ICON_FAST_REWIND ("transport_ctrl_fast-rewind.png")
#define GUI_VIEWCUBE_ROTATION_ARROW ("viewcube_rotation_arrow.png")
#define GUI_VIEWCUBE_UP_ARROW ("viewcube_up_arrow.png")


namespace megamol {
namespace gui {

    /********** Additional Global ImGui Operators ****************************/

    namespace {

        bool operator==(const ImVec2& left, const ImVec2& right) {
            return ((left.x == right.x) && (left.y == right.y));
        }

        bool operator!=(const ImVec2& left, const ImVec2& right) {
            return !(left == right);
        }

    } // namespace


    /********** Global Unique ID *********************************************/

    /// ! Do not directly change
    extern ImGuiID gui_generated_uid;

    inline ImGuiID GenerateUniqueID(void) {
        return (++gui_generated_uid);
    }


    /********** Global ImGui Context Pointer Counter *************************/

    // Only accessed by possible multiple instances of GUIManager
    extern unsigned int gui_context_count;


    /********** Global Resource Paths ****************************************/

    // Resource paths set by GUIManager
    extern std::vector<std::string> gui_resource_paths;


    /********** Global GUI Scaling Factor ************************************/

    // Forward declaration
    class GUIManager;

    class GUIScaling {
    public:
        friend class GUIManager;

        GUIScaling() = default;
        ~GUIScaling() = default;

        float Get(void) const {
            return this->scale;
        }

        float TransitionFactor(void) const {
            return (this->scale / this->last_scale);
        }

        bool PendingChange(void) const {
            return this->pending_change;
        }

    private:
        bool ConsumePendingChange(void) {
            bool current_pending_change = this->pending_change;
            this->pending_change = false;
            return current_pending_change;
        }

        void Set(float s) {
            this->last_scale = this->scale;
            this->scale = std::max(0.0f, s);
            if (this->scale != this->last_scale) {
                this->pending_change = true;
            }
        }

        bool pending_change = false;
        float scale = 1.0f;
        float last_scale = 1.0f;
    };

    extern GUIScaling gui_scaling;


    /********** Types ********************************************************/

    // Forward declaration
    class CallSlot;
    class InterfaceSlot;
    typedef std::shared_ptr<megamol::gui::CallSlot> CallSlotPtr_t;
    typedef std::shared_ptr<megamol::gui::InterfaceSlot> InterfaceSlotPtr_t;

    /** Available ImGui APIs */
    enum class GUIImGuiAPI { NONE, OPEN_GL };

    /** Hotkey Data Types (exclusively for configurator) */
    enum HotkeyIndex : size_t {
        MODULE_SEARCH = 0,
        PARAMETER_SEARCH = 1,
        DELETE_GRAPH_ITEM = 2,
        SAVE_PROJECT = 3,
        INDEX_COUNT = 4
    };

    struct HotkeyData_t {
        megamol::core::view::KeyCode keycode;
        bool is_pressed = false;
    };

    typedef std::array<megamol::gui::HotkeyData_t, megamol::gui::HotkeyIndex::INDEX_COUNT> HotkeyArray_t;

    typedef megamol::core::param::AbstractParamPresentation::Presentation Present_t;
    typedef megamol::core::param::AbstractParamPresentation::ParamType ParamType_t;
    typedef std::map<int, std::string> EnumStorage_t;

    typedef std::array<float, 5> FontScalingArray_t;

    /* Data type holding a pair of uids. */
    typedef std::vector<ImGuiID> UIDVector_t;
    typedef std::pair<ImGuiID, ImGuiID> UIDPair_t;
    typedef std::vector<UIDPair_t> UIDPairVector_t;

    typedef std::pair<std::string, std::string> StrPair_t;
    typedef std::vector<StrPair_t> StrPairVector_t;

    /* Data type holding current group uid and group name pairs. */
    typedef std::pair<ImGuiID, std::string> GraphGroupPair_t;
    typedef std::vector<megamol::gui::GraphGroupPair_t> GraphGroupPairVector_t;

    enum PresentPhase : size_t { INTERACTION = 0, RENDERING = 1 };

    /* Data type holding information of graph canvas. */
    typedef struct _canvas_ {
        ImVec2 position;      // in
        ImVec2 size;          // in
        ImVec2 scrolling;     // in
        float zooming;        // in
        ImVec2 offset;        // in
        ImFont* gui_font_ptr; // in
    } GraphCanvas_t;

    /* Data type holding information on graph item interaction. */
    typedef struct _interact_state_ {

        ImGuiID button_active_uid;  // in out
        ImGuiID button_hovered_uid; // in out
        bool process_deletion;      // out

        ImGuiID group_selected_uid; // in out
        ImGuiID group_hovered_uid;  // in out
        bool group_layout;          // out

        UIDVector_t modules_selected_uids;               // in out
        ImGuiID module_hovered_uid;                      // in out
        UIDPairVector_t modules_add_group_uids;          // out
        UIDVector_t modules_remove_group_uids;           // out
        bool modules_layout;                             // out
        StrPairVector_t module_rename;                   // out
        vislib::math::Ternary module_graphentry_changed; // out
        ImVec2 module_param_child_position;              // out
        bool module_show_label;                          // in

        ImGuiID call_selected_uid;  // in out
        ImGuiID call_hovered_uid;   // in out
        bool call_show_label;       // in
        bool call_show_slots_label; // in

        ImGuiID slot_dropped_uid; // in out

        ImGuiID callslot_selected_uid;       // in out
        ImGuiID callslot_hovered_uid;        // in out
        UIDPair_t callslot_add_group_uid;    // in out
        UIDPair_t callslot_remove_group_uid; // in out
        CallSlotPtr_t callslot_compat_ptr;   // in
        bool callslot_show_label;            // in

        ImGuiID interfaceslot_selected_uid;          // in out
        ImGuiID interfaceslot_hovered_uid;           // in out
        InterfaceSlotPtr_t interfaceslot_compat_ptr; // in

        bool parameters_extended_mode; // in

        bool graph_is_running; // in

    } GraphItemsInteract_t;

    /* Data type holding shared state of graph items. */
    typedef struct _graph_item_state_ {
        megamol::gui::GraphCanvas_t canvas;          // (see above)
        megamol::gui::GraphItemsInteract_t interact; // (see above)
        megamol::gui::HotkeyArray_t hotkeys;         // in out
        megamol::gui::GraphGroupPairVector_t groups; // in
    } GraphItemsState_t;

    /* Data type holding shared state of graphs. */
    typedef struct _graph_state_ {
        FontScalingArray_t graph_zoom_font_scalings; // in
        float graph_width;                           // in
        bool show_parameter_sidebar;                 // in
        megamol::gui::HotkeyArray_t hotkeys;         // in out
        ImGuiID graph_selected_uid;                  // out
        bool graph_delete;                           // out
        bool configurator_graph_save;                // out
        bool global_graph_save;                      // out
        ImGuiID new_running_graph_uid;               // out
    } GraphState_t;

    enum class HeaderType { MODULE_GROUP, MODULE, PARAMETER_GROUP };


    /********** GUIUtils *****************************************************/

    /**
     * Static GUI utility functions.
     */
    class GUIUtils {
    public:
        /** Extract string enclosed in predefined tags. */
        static std::string ExtractTaggedString(
            const std::string& str, const std::string& start_tag, const std::string& end_tag) {
            std::string return_str;
            auto start_idx = str.find(start_tag);
            if (start_idx != std::string::npos) {
                auto end_idx = str.find(end_tag, start_idx);
                if ((end_idx != std::string::npos) && (start_idx < end_idx)) {
                    start_idx += std::string(start_tag).length();
                    return_str = str.substr(start_idx, (end_idx - start_idx));
                }
            }
            return return_str;
        }

        /** Decode string from UTF-8. */
        static bool Utf8Decode(std::string& str) {
            vislib::StringA dec_tmp;
            if (vislib::UTF8Encoder::Decode(dec_tmp, vislib::StringA(str.c_str()))) {
                str = std::string(dec_tmp.PeekBuffer());
                return true;
            }
            return false;
        }

        /** Encode string into UTF-8. */
        static bool Utf8Encode(std::string& str) {
            vislib::StringA dec_tmp;
            if (vislib::UTF8Encoder::Encode(dec_tmp, vislib::StringA(str.c_str()))) {
                str = std::string(dec_tmp.PeekBuffer());
                return true;
            }
            return false;
        }

        /**
         * Enable/Disable read only widget style.
         */
        static void ReadOnlyWigetStyle(bool set) {
            if (set) {
                ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
                ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
            } else {
                ImGui::PopItemFlag();
                ImGui::PopStyleVar();
            }
        }

        /**
         * Returns true if search string is found in source as a case insensitive substring.
         *
         * @param source   The string to search in.
         * @param search   The string to search for in the source.
         */
        static bool FindCaseInsensitiveSubstring(const std::string& source, const std::string& search) {
            if (search.empty())
                return true;
            auto it = std::search(source.begin(), source.end(), search.begin(), search.end(),
                [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); });
            return (it != source.end());
        }

        /**
         * Returns true if both strings equal each other case insensitively.
         *
         * @param source   One string.
         * @param search   Second string.
         */
        static bool CaseInsensitiveStringCompare(std::string const& str1, std::string const& str2) {
            return ((str1.size() == str2.size()) &&
                    std::equal(str1.begin(), str1.end(), str2.begin(), [](char const& c1, char const& c2) {
                        return (c1 == c2 || std::toupper(c1) == std::toupper(c2));
                    }));
        }

        /*
         * Draw collapsing group header.
         */
        static bool GroupHeader(megamol::gui::HeaderType type, const std::string& name, std::string& inout_search,
            ImGuiID override_header_state = GUI_INVALID_ID) {
            if (ImGui::GetCurrentContext() == nullptr) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "[GUI] No ImGui context available. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
                return false;
            }

            // Determine header state and change color depending on active parameter search
            auto headerId = ImGui::GetID(name.c_str());
            auto headerState = override_header_state;
            if (headerState == GUI_INVALID_ID) {
                headerState = ImGui::GetStateStorage()->GetInt(headerId, 0); // 0=close 1=open
            }

            int pop_style_color_number = 0;
            if (type == megamol::gui::HeaderType::MODULE_GROUP) {
                ImGui::PushStyleColor(ImGuiCol_Header, GUI_COLOR_GROUP_HEADER);
                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, GUI_COLOR_GROUP_HEADER_HIGHLIGHT);
                auto header_active_color = GUI_COLOR_GROUP_HEADER_HIGHLIGHT;
                header_active_color.w *= 0.75f;
                ImGui::PushStyleColor(ImGuiCol_HeaderActive, header_active_color);
                pop_style_color_number += 3;
            }

            bool searched = true;
            if (!inout_search.empty()) {
                headerState = 1;
                searched = GUIUtils::FindCaseInsensitiveSubstring(name, inout_search);
                if (!searched) {
                    auto header_col = ImGui::GetStyleColorVec4(ImGuiCol_Header);
                    header_col.w *= 0.25;
                    ImGui::PushStyleColor(ImGuiCol_Header, header_col);
                    auto header_hovered_col = ImGui::GetStyleColorVec4(ImGuiCol_HeaderHovered);
                    header_hovered_col.w *= 0.25;
                    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, header_hovered_col);
                    auto text_col = ImGui::GetStyleColorVec4(ImGuiCol_Text);
                    text_col.w *= 0.25;
                    ImGui::PushStyleColor(ImGuiCol_Text, text_col);
                    pop_style_color_number += 3;
                } else {
                    // Show all below when given name is part of the search
                    inout_search.clear();
                }
            }
            ImGui::GetStateStorage()->SetInt(headerId, headerState);
            bool header_open = ImGui::CollapsingHeader(name.c_str(), nullptr);
            ImGui::PopStyleColor(pop_style_color_number);

            // Keep following elements open for one more frame to propagate override changes to headers below.
            if (override_header_state == 0) {
                header_open = true;
            }
            return header_open;
        }

    private:
        GUIUtils(void) = default;
        ~GUIUtils(void) = default;
    };

} // namespace gui
} // namespace megamol

#endif // MEGAMOL_GUI_GUIUTILS_INCLUDED
