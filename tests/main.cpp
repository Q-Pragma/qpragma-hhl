/**
 * Collect all tests and run tests
 * using google tests
 */

// Include Google tests
#include "gtest/gtest.h"


// Run tests
int main(int argc, char ** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
