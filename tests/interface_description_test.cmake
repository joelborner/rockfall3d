# Exercise the executable, independently of any plugin or optional model fixture.
file(MAKE_DIRECTORY "${WORK_DIR}")
file(GLOB_RECURSE before LIST_DIRECTORIES true "${WORK_DIR}/*")
if(before)
    message(FATAL_ERROR "Interface command test requires an empty working directory: ${WORK_DIR}")
endif()

execute_process(
    COMMAND "${EXECUTABLE}" --describe-interface
    WORKING_DIRECTORY "${WORK_DIR}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE description
    ERROR_VARIABLE diagnostics
    TIMEOUT 15
)
if(NOT result STREQUAL "0" OR NOT diagnostics STREQUAL "")
    message(FATAL_ERROR "Interface command failed (${result}): ${diagnostics}")
endif()

# JSON parsing rejects stray logs. The generated catalogue is authoritative;
# no particular optional model is required or forbidden by this test.
string(JSON schema_type TYPE "${description}" interface_schema_version)
string(JSON schema_version GET "${description}" interface_schema_version)
if(NOT schema_type STREQUAL "NUMBER" OR NOT schema_version STREQUAL "1")
    message(FATAL_ERROR "Unexpected interface schema: ${schema_version}")
endif()
string(JSON engine_id GET "${description}" engine id)
if(NOT engine_id STREQUAL "rockfall3d")
    message(FATAL_ERROR "Unexpected engine ID: ${engine_id}")
endif()
foreach(key version build_id description_id)
    string(JSON value_type TYPE "${description}" engine "${key}")
    string(JSON value GET "${description}" engine "${key}")
    if(NOT value_type STREQUAL "STRING" OR value STREQUAL "")
        message(FATAL_ERROR "Missing engine identity field: ${key}")
    endif()
endforeach()
foreach(key managed_domains sections)
    string(JSON value_type TYPE "${description}" "${key}")
    if(NOT value_type STREQUAL "ARRAY")
        message(FATAL_ERROR "Expected an array: ${key}")
    endif()
endforeach()

file(READ "${EXPECTED_DESCRIPTION}" expected)
string(REPLACE "\r\n" "\n" expected "${expected}")
string(REPLACE "\r\n" "\n" description "${description}")
if(NOT description STREQUAL expected)
    message(FATAL_ERROR "Executable stdout does not exactly match its generated catalogue")
endif()

foreach(arguments IN ITEMS "" "--describe-interface;unexpected")
    execute_process(
        COMMAND "${EXECUTABLE}" ${arguments}
        WORKING_DIRECTORY "${WORK_DIR}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE diagnostics
        TIMEOUT 15
    )
    if(result STREQUAL "0" OR NOT output STREQUAL "" OR NOT diagnostics MATCHES "usage:")
        message(FATAL_ERROR "Invalid argument handling failed for '${arguments}': ${result}, ${diagnostics}")
    endif()
endforeach()

file(GLOB_RECURSE after LIST_DIRECTORIES true "${WORK_DIR}/*")
if(after)
    message(FATAL_ERROR "Interface/argument commands created files: ${after}")
endif()
