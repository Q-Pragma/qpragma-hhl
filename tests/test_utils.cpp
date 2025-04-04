/**
 * This file ensure the state preparation works as
 * expected and produce the right state
 */

// Include headers
#include "gtest/gtest.h"

#include "qpragma.h"
#include "qpragma/hhl.hpp"


// Ensure normalization works
TEST(utils, normalize) {
    /*
     * Test 1
     */

    // Create array
    std::array<double, 2UL> arr_1 { 1., 1. };
    qpragma::hhl::utils::normalize<1UL>(arr_1);

    // Check result
    ASSERT_DOUBLE_EQ(arr_1[0UL], 1. / std::sqrt(2.));
    ASSERT_DOUBLE_EQ(arr_1[1UL], 1. / std::sqrt(2.));

    /*
     * Test 2
     */

    // Create array
    std::array<double, 4UL> arr_2 { 1., 0., 0., 0. };
    qpragma::hhl::utils::normalize<2UL>(arr_2);

    // Check result
    ASSERT_DOUBLE_EQ(arr_2[0UL], 1.);
    ASSERT_DOUBLE_EQ(arr_2[1UL], 0.);
    ASSERT_DOUBLE_EQ(arr_2[2UL], 0.);
    ASSERT_DOUBLE_EQ(arr_2[3UL], 0.);

}
