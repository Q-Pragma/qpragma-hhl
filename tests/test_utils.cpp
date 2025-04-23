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

// Ensure sign function works properly
TEST(utils, sign) {
    /* Test on negative values */
    double val = -0.05;
    ASSERT_DOUBLE_EQ(qpragma::hhl::utils::sign(val), -1.);
    val = -4.;
    ASSERT_DOUBLE_EQ(qpragma::hhl::utils::sign(val), -1.);

    /* Test on positive values */
    val = 0.05;
    ASSERT_DOUBLE_EQ(qpragma::hhl::utils::sign(val), 1.);
    val = 4.;
    ASSERT_DOUBLE_EQ(qpragma::hhl::utils::sign(val), 1.);
}

// Ensure get_means function works properly
TEST(utils, get_means) {
    // /!\ Bits are read in reverse order

    /* [0=0b00, 1=0b01] -- > [1/2, 0] */
    {
        constexpr uint64_t SIZE = 2UL;
        std::vector<uint64_t> eigenvals {0, 1};
        std::array<double, SIZE> means = qpragma::hhl::utils::get_means<SIZE>(eigenvals);
        std::array<double, SIZE> expected {1./2., 0.};
        for (uint64_t idx = 0UL ; idx < SIZE ; ++idx) {
            ASSERT_DOUBLE_EQ(means[idx], expected[idx]);
        }
    }

    /* [0=0b000, 4=0b100, 5=0b101] -- > [1/3, 0, 2/3] */
    {
        constexpr uint64_t SIZE = 3UL;
        std::vector<uint64_t> eigenvals {0, 4, 5};
        std::array<double, SIZE> means = qpragma::hhl::utils::get_means<SIZE>(eigenvals);
        std::array<double, SIZE> expected {1./3., 0. , 2./3.};
        for (uint64_t idx = 0UL ; idx < SIZE ; ++idx) {
            ASSERT_DOUBLE_EQ(means[idx], expected[idx]);
        }
    }

    /* [12=0b1100, 4=0b0100, 6=0b0110] -- > [0, 1/3, 1, 1/3] */
    {
        constexpr uint64_t SIZE = 4UL;
        std::vector<uint64_t> eigenvals {12, 4, 6};
        std::array<double, SIZE> means = qpragma::hhl::utils::get_means<SIZE>(eigenvals);
        std::array<double, SIZE> expected {0., 1./3., 1., 1./3.};
        for (uint64_t idx = 0UL ; idx < SIZE ; ++idx) {
            ASSERT_DOUBLE_EQ(means[idx], expected[idx]);
        }
    }

}

// Ensure bin_to_double function works properly
TEST(utils, bin_to_double) {
    /* 0 = 0b00 -- > 0. */
    {
        constexpr uint64_t SIZE = 2UL;
        uint64_t val = 0;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), 0.);
    }

    /* 0 = 0b01 -- > 0.5 */
    {
        constexpr uint64_t SIZE = 2UL;
        uint64_t val = 1;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), 0.5);
    }

    /* 0 = 0b10 -- > 1. */
    {
        constexpr uint64_t SIZE = 2UL;
        uint64_t val = 2;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), 1);
    }

    /* 0 = 0b11 -- > -0.5 */
    {
        constexpr uint64_t SIZE = 2UL;
        uint64_t val = 3;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), -0.5);
    }

    /* 1 = 0b001 -- > 0.25 */
    {
        constexpr uint64_t SIZE = 3UL;
        uint64_t val = 1;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), 0.25);
    }

    /* 5 = 0b101 -- > -0.75 */
    {
        constexpr uint64_t SIZE = 3UL;
        uint64_t val = 5;
        ASSERT_DOUBLE_EQ(qpragma::hhl::utils::bin_to_double(SIZE, val), -0.75);
    }
}
