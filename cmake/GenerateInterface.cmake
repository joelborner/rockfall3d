set(ENGINE_BUILD_ID "source-archive")
set(ENGINE_SOURCE_REVISION "null")
set(ENGINE_SOURCE_DIRTY "null")

if(GIT_EXECUTABLE AND EXISTS "${GIT_EXECUTABLE}")
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${ROCKFALL3D_SOURCE_DIR}" rev-parse HEAD
        RESULT_VARIABLE revision_result
        OUTPUT_VARIABLE revision
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(revision_result EQUAL 0 AND revision MATCHES "^[0-9a-f]+$")
        set(ENGINE_BUILD_ID "${revision}")
        set(ENGINE_SOURCE_REVISION "\"${revision}\"")
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${ROCKFALL3D_SOURCE_DIR}"
                status --porcelain --untracked-files=normal
            RESULT_VARIABLE status_result
            OUTPUT_VARIABLE source_status
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(status_result EQUAL 0)
            set(ENGINE_SOURCE_DIRTY "false")
            if(NOT source_status STREQUAL "")
                set(ENGINE_SOURCE_DIRTY "true")
                string(APPEND ENGINE_BUILD_ID "-dirty")
            endif()
        endif()
    endif()
endif()

set(interface_source "${ROCKFALL3D_SOURCE_DIR}/interface/description.json.in")
file(SHA256 "${interface_source}" ENGINE_DESCRIPTION_HASH)
file(MAKE_DIRECTORY "${ROCKFALL3D_GENERATED_DIR}")
configure_file("${interface_source}"
    "${ROCKFALL3D_GENERATED_DIR}/interface-description.json" @ONLY)
file(READ "${ROCKFALL3D_GENERATED_DIR}/interface-description.json" INTERFACE_JSON)
string(JSON schema_version GET "${INTERFACE_JSON}" interface_schema_version)
if(NOT schema_version EQUAL 1)
    message(FATAL_ERROR "Unsupported interface_schema_version: ${schema_version}")
endif()
configure_file("${ROCKFALL3D_SOURCE_DIR}/cmake/interface_description.hpp.in"
    "${ROCKFALL3D_GENERATED_DIR}/interface_description.hpp" @ONLY)
