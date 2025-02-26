# Define C++ version used
find_program(CXX clang++)
set(CMAKE_CXX_COMPILER ${CXX})
set(CMAKE_CXX_STANDARD 20)

# Set environment variables
set(PROJECT_DIR ${CMAKE_SOURCE_DIR})
set(SRC_DIR ${PROJECT_DIR}/src)
set(INCLUDE_DIR ${PROJECT_DIR}/include)
set(BUILD_DIR ${CMAKE_BINARY_DIR})
set(TESTS_DIR ${PROJECT_DIR}/tests)
