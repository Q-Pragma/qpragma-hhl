# ####################### #
# Configure cmake project #
# ####################### #

# Define project
project(qpragma-hhl LANGUAGES CXX)
cmake_minimum_required(VERSION 3.10.2)

# Define useful variables
include("cmake/environ.cmake")


# #################### #
# Q-Pragma HHL Library #
# #################### #

# Define C++ files
include_directories(${INCLUDE_DIR})

set(qpragma-hhl-cpp
        ${SRC_DIR}/utils.cpp)

set(qpragma-hhl-headers
        ${INCLUDE_DIR}/qpragma/hhl.hpp
        ${INCLUDE_DIR}/qpragma/hhl/core.hpp
        ${INCLUDE_DIR}/qpragma/hhl/observables.hpp
        ${INCLUDE_DIR}/qpragma/hhl/observables.ipp
        ${INCLUDE_DIR}/qpragma/hhl/simulation.hpp
        ${INCLUDE_DIR}/qpragma/hhl/stateprep.hpp
        ${INCLUDE_DIR}/qpragma/hhl/stateprep.ipp
        ${INCLUDE_DIR}/qpragma/hhl/utils.hpp
        ${INCLUDE_DIR}/qpragma/hhl/utils.ipp)

# Define compilation rules
add_library(qpragma-hhl SHARED ${qpragma-hhl-cpp})
target_link_libraries(qpragma-hhl qpragma)
set_target_properties(
    qpragma-hhl
        PROPERTIES
            PUBLIC_HEADER "${qpragma-hhl-headers}"
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS -fplugin=qpragma-plugin.so
)


# ####################### #
# Q-Pragma HHL executable #
# ####################### #

# Define compilation rules
add_executable(qpragma-hhl-exec ${SRC_DIR}/main.cpp)
target_link_libraries(qpragma-hhl-exec qpragma-hhl qpragma qpragma-newlinalg qatnewlinalg boost_program_options)
set_target_properties(
    qpragma-hhl-exec
    PROPERTIES
        LINKER_LANGUAGE CXX
        COMPILE_FLAGS -fplugin=qpragma-plugin.so
        RUNTIME_OUTPUT_NAME qpragma-hhl
)


# ################################## #
# Install HHL library and executable #
# ################################## #

# Install library
install(TARGETS qpragma-hhl
        LIBRARY DESTINATION usr/lib64)
        # PUBLIC_HEADER DESTINATION usr/include)

# Install executable
install(TARGETS qpragma-hhl-exec
        RUNTIME DESTINATION usr/bin)


# ######################### #
# Run C++ tests using gtest #
# ######################### #

add_subdirectory(tests)
