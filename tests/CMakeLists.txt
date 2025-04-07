# ################################ #
# Test Q-Pragma HHL implementation #
# ################################ #

# List C++ files
set(cpp-test-files
        ${TESTS_DIR}/main.cpp
        ${TESTS_DIR}/test_hhl.cpp
        ${TESTS_DIR}/test_utils.cpp)


# Use Q-Pragma plugin
add_executable(qpragma-hhl-tests EXCLUDE_FROM_ALL ${cpp-test-files})
target_link_libraries(qpragma-hhl-tests qpragma qpragma-newlinalg qatnewlinalg gtest)
set_target_properties(
    qpragma-hhl-tests
    PROPERTIES
        LINKER_LANGUAGE CXX
        COMPILE_FLAGS -fplugin=qpragma-plugin.so
)


# Add C++ tests target
add_custom_target(cpp_tests
    DEPENDS qpragma-hhl-tests
    COMMAND $<TARGET_FILE:qpragma-hhl-tests>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)


# Create check target
add_custom_target(check DEPENDS cpp_tests)
