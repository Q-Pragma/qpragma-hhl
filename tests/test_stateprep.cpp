// Unitary test using gtest
#include <string>
#include "gtest/gtest.h"
#include "qpragma.h"
#include "qpragma/hhl.hpp"

using namespace qpragma;

constexpr uint64_t NB_SHOTS = 1000UL;
// Rough confidence interval
constexpr double ERROR_RATIO = 0.05;


/* Tests for the state preparation functions */

// Tests for KPTree (only in probability with this form of QPU)

TEST(KPTree, 1qubit) {
    constexpr uint64_t SIZE = 1UL;

    quint_t<SIZE> qreg;

    std::vector<std::array<double, (1 << SIZE)>> init_arrays {
        {1., 0.},
        {0., 1.},
        {1., 1.},
        {0.34, 0.91},
        {0.99, -0.01}
    };

    for (auto init_array : init_arrays) {
        hhl::utils::normalize<SIZE>(init_array);
        std::array<double, (1 << SIZE)> result;
        result.fill(0.);
        for (uint64_t shot = 0UL ; shot < NB_SHOTS ; ++shot) {
            (hhl::stateprep::kp_tree<SIZE>(init_array))(qreg);
            ++result[measure_and_reset(qreg)];
        }

        for (uint64_t idx = 0UL ; idx < (1 << SIZE) ; ++idx) {
            double expected_val = init_array[idx] * init_array[idx];
            ASSERT_NEAR(result[idx] / NB_SHOTS, expected_val, ERROR_RATIO);
        }
    }
}

TEST(KPTree, 2qubits) {
    constexpr uint64_t SIZE = 2UL;

    quint_t<SIZE> qreg;

    std::vector<std::array<double, (1 << SIZE)>> init_arrays {
        {1., 0., 0., 0.},
        {0., 1., 0., 0.},
        {0., 0., 1., 0.},
        {0., 0., 0., 1.},
        {1., 1., 1., 1.},
        {0.34, 0.91, 0.23, -0.1},
        {0.99, 0.01, 0.2, 0.},
        {1., 2., 3., 4.}
    };

    for (auto init_array : init_arrays) {
        hhl::utils::normalize<SIZE>(init_array);
        std::array<double, (1 << SIZE)> result;
        result.fill(0.);
        for (uint64_t shot = 0UL ; shot < NB_SHOTS ; ++shot) {
            (hhl::stateprep::kp_tree<SIZE>(init_array))(qreg);
            ++result[measure_and_reset(qreg)];
        }

        for (uint64_t idx = 0UL ; idx < (1 << SIZE) ; ++idx) {
            double expected_val = init_array[idx] * init_array[idx];
            ASSERT_NEAR(result[idx] / NB_SHOTS, expected_val, ERROR_RATIO);
        }
    }
}

TEST(KPTree, 3qubits) {
    constexpr uint64_t SIZE = 3UL;

    quint_t<SIZE> qreg;

    std::vector<std::array<double, (1 << SIZE)>> init_arrays {
        {1., 0., 0., 0., 0., 0., 0., 0.},
        {0., 1., 0., 0., 0., 0., 0., 0.},
        {0., 0., 1., 0., 0., 0., 0., 0.},
        {0., 0., 0., 1., 0., 0., 0., 0.},
        {0., 0., 0., 0., 1., 0., 0., 0.},
        {0., 0., 0., 0., 0., 1., 0., 0.},
        {0., 0., 0., 0., 0., 0., 1., 0.},
        {0., 0., 0., 0., 0., 0., 0., 1.},
        {1., 1., 1., 1., 1., 1., 1., 1.},
        {0.34, 0.91, 0.23, -0.1, 0., 0., 0.,0.},
        {0.99, 0.01, 0.2, 0., 0.2, 0.34, 0.1, 0.88},
        {1., 2., 3., 4., 5., 6., 7., 8.},
        {3., 1.4, 1.5, 9.2, 6.5, 3.5, 8.9, 7.9}
    };

    for (auto init_array : init_arrays) {
        hhl::utils::normalize<SIZE>(init_array);
        std::array<double, (1 << SIZE)> result;
        result.fill(0.);
        for (uint64_t shot = 0UL ; shot < NB_SHOTS ; ++shot) {
            (hhl::stateprep::kp_tree<SIZE>(init_array))(qreg);
            ++result[measure_and_reset(qreg)];
        }

        for (uint64_t idx = 0UL ; idx < (1 << SIZE) ; ++idx) {
            double expected_val = init_array[idx] * init_array[idx];
            ASSERT_NEAR(result[idx] / NB_SHOTS, expected_val, ERROR_RATIO);
        }
    }
}