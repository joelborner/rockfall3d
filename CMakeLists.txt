cmake_minimum_required(VERSION 3.21)

project(rockfall3d LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ROCKFALL3D_WITH_GDAL "Build with GDAL-based raster and shapefile support" ON)

add_executable(rockfall3d
    rockfall3d/main.cpp
    rockfall3d/rock.cpp
    rockfall3d/terrain.cpp
    rockfall3d/roughness.cpp
    rockfall3d/time_stepper.cpp
    rockfall3d/config_io.hpp
    rockfall3d/contact_wmat.hpp
    rockfall3d/convex_polyhedron.cpp
    rockfall3d/convex_polyhedron.hpp
    rockfall3d/deposition_output.hpp
    rockfall3d/fast_io.hpp
    rockfall3d/fragmentation.cpp
    rockfall3d/fragmentation.hpp
    rockfall3d/jor_solver.hpp
    rockfall3d/math_alias.hpp
    rockfall3d/polygon_utils.hpp
    rockfall3d/raster_output.hpp
    rockfall3d/rock.hpp
    rockfall3d/rotational_integrator.hpp
    rockfall3d/rolling_friction.hpp
    rockfall3d/roughness.hpp
    rockfall3d/sor_solver.hpp
    rockfall3d/terrain.hpp
    rockfall3d/terrain_export.hpp
    rockfall3d/terrain_interactions.hpp
    rockfall3d/terrain_io.hpp
    rockfall3d/time_stepper.hpp
    rockfall3d/trajectory_format.hpp
    rockfall3d/trajectory_recorder.hpp
)

target_include_directories(rockfall3d
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
        ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d/external/quickhull/include
)

find_package(Eigen3 3.4 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
    target_link_libraries(rockfall3d PRIVATE Eigen3::Eigen)
else()
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Dense
        DOC "Path to the Eigen headers"
    )
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR "Eigen3 was not found. Install Eigen or configure EIGEN3_INCLUDE_DIR.")
    endif()

    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}"
    )
    target_link_libraries(rockfall3d PRIVATE Eigen3::Eigen)
endif()

if(BUILD_TESTING)
    add_executable(rotational_integrator_tests
        tests/rotational_integrator_tests.cpp
        rockfall3d/math_alias.hpp
        rockfall3d/rotational_integrator.hpp
    )
    target_include_directories(rotational_integrator_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
    )
    target_link_libraries(rotational_integrator_tests PRIVATE Eigen3::Eigen)
    add_test(NAME rotational_integrator_tests COMMAND rotational_integrator_tests)

    add_executable(rolling_friction_tests
        tests/rolling_friction_tests.cpp
        rockfall3d/jor_solver.hpp
        rockfall3d/math_alias.hpp
        rockfall3d/rolling_friction.hpp
        rockfall3d/sor_solver.hpp
    )
    target_include_directories(rolling_friction_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
    )
    target_link_libraries(rolling_friction_tests PRIVATE Eigen3::Eigen)
    add_test(NAME rolling_friction_tests COMMAND rolling_friction_tests)

    add_executable(contact_solver_tests tests/contact_solver_tests.cpp)
    target_include_directories(contact_solver_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
    )
    target_link_libraries(contact_solver_tests PRIVATE Eigen3::Eigen)
    add_test(NAME contact_solver_tests COMMAND contact_solver_tests)

    add_executable(terrain_interaction_tests
        tests/terrain_interaction_tests.cpp
        rockfall3d/math_alias.hpp
        rockfall3d/terrain_interactions.hpp
    )
    target_include_directories(terrain_interaction_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
    )
    target_link_libraries(terrain_interaction_tests PRIVATE Eigen3::Eigen)
    add_test(NAME terrain_interaction_tests COMMAND terrain_interaction_tests)

    add_executable(fragmentation_tests
        tests/fragmentation_tests.cpp
        rockfall3d/convex_polyhedron.cpp
        rockfall3d/fragmentation.cpp
        rockfall3d/rock.cpp
    )
    target_include_directories(fragmentation_tests
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d
            ${CMAKE_CURRENT_SOURCE_DIR}/rockfall3d/external/quickhull/include
    )
    target_link_libraries(fragmentation_tests PRIVATE Eigen3::Eigen)
    add_test(NAME fragmentation_tests COMMAND fragmentation_tests)
endif()

if(ROCKFALL3D_WITH_GDAL)
    find_package(GDAL REQUIRED)
    target_compile_definitions(rockfall3d PRIVATE WITH_GDAL)

    if(TARGET GDAL::GDAL)
        target_link_libraries(rockfall3d PRIVATE GDAL::GDAL)
    else()
        target_include_directories(rockfall3d PRIVATE ${GDAL_INCLUDE_DIRS})
        target_link_libraries(rockfall3d PRIVATE ${GDAL_LIBRARIES})
    endif()
endif()

if(MSVC)
    target_compile_options(rockfall3d PRIVATE /W3 /permissive- /EHsc)
else()
    target_compile_options(rockfall3d PRIVATE -Wall -Wextra -Wpedantic)
endif()
