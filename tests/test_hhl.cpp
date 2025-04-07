// Unitary test using gtest
#include <string>
#include "gtest/gtest.h"
#include "qpragma.h"
#include "qpragma/hhl.hpp"

using namespace qpragma;


double ERROR_RATIO = 0.1;

/* 1 qubit - 2x2 matrices */


TEST(HHL, m2x2_I) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 I */
    hhl::observables::PauliTerm<SIZE> term1 { "I", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {1., 0.};

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
    {
        // b = |1>
        std::array<double, 1UL << SIZE> init_array {0., 1.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 1.};
    
        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);
    
        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m2x2_X) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 X */
    hhl::observables::PauliTerm<SIZE> term1 { "X", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 1.};

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
    {
        // b = |1>
        std::array<double, 1UL << SIZE> init_array {0., 1.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {1., 0.};
    
        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);
    
        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m2x2_Y) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 Y */
    hhl::observables::PauliTerm<SIZE> term1 { "Y", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 1.};

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
    {
        // b = |1>
        std::array<double, 1UL << SIZE> init_array {0., 1.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {1., 0.};
    
        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);
    
        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m2x2_Z) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 Z */
    hhl::observables::PauliTerm<SIZE> term1 { "Z", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {1., 0.};

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
    {
        // b = |1>
        std::array<double, 1UL << SIZE> init_array {0., 1.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 1.};
    
        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);
    
        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m2x2_usecase1) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 I + 0.25 X */
    hhl::observables::PauliTerm<SIZE> term1 { "I", 0.5 };
    hhl::observables::PauliTerm<SIZE> term2 { "X", 0.25};
    std::vector term_vect { term1, term2 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0.8, 0.2};

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
    {
        // b = |1>
        std::array<double, 1UL << SIZE> init_array {0., 1.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0.2, 0.8};
    
        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);
    
        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}