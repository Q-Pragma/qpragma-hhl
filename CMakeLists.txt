# ####################### #
# Configure cmake project #
# ####################### #

# Define project
project(qpragma-hhl LANGUAGES CXX)
cmake_minimum_required(VERSION 3.10.2)

# Define useful variables
include("cmake/environ.cmake")


# ######################## #
# Q-Pragma shor executable #
# ######################## #

# Define C++ files
include_directories(${INCLUDE_DIR})

set(qpragma-hhl-cpp
        ${SRC_DIR}/observable.cpp)

set(qpragma-hhl-headers
        # ${INCLUDE_DIR}/qpragma/hhl.hpp
        ${INCLUDE_DIR}/qpragma/hhl/observables.hpp
        ${INCLUDE_DIR}/qpragma/hhl/simulations.hpp
        ${INCLUDE_DIR}/qpragma/hhl/stateprep.hpp)

# Define compilation rules
add_executable(qpragma-hhl ${qpragma-hhl-cpp} ${SRC_DIR}/main.cpp)
target_link_libraries(qpragma-hhl qpragma qpragma-newlinalg qatnewlinalg boost_program_options)
set_target_properties(
    qpragma-hhl PROPERTIES PUBLIC_HEADER "${qpragma-hhl-headers}"
                            LINKER_LANGUAGE CXX
                            COMPILE_FLAGS -fplugin=qpragma-plugin.so)

# Install
install(TARGETS qpragma-hhl
        RUNTIME DESTINATION usr/bin)
