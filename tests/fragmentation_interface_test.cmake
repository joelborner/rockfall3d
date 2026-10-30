file(READ "${DESCRIPTION}" description)
string(JSON section_id GET "${description}" sections 0 id)
if(NOT section_id STREQUAL "fragmentation")
    message(FATAL_ERROR "Fragmentation-enabled build must advertise fragmentation")
endif()

string(JSON element GET "${description}" sections 0 xml element)
set(attributes "")
set(values "")
set(defaults "")
foreach(index RANGE 0 3)
    string(JSON attribute GET "${description}" sections 0 parameters ${index} xml attribute)
    string(JSON initial GET "${description}" sections 0 parameters ${index} initial_value)
    string(APPEND attributes " ${attribute}=\"${initial}\"")
    list(APPEND values "${initial}")
    if(index GREATER 1)
        string(JSON default GET "${description}" sections 0 parameters ${index} parser_default)
        list(APPEND defaults "${default}")
    endif()
endforeach()
execute_process(
    COMMAND "${PARSER_TEST}" "${WORK_DIR}" "<${element}${attributes}/>" ${values} ${defaults}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE diagnostics
    TIMEOUT 30
)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Fragmentation catalogue/parser mismatch: ${output}${diagnostics}")
endif()
