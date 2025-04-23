// Unitary test using gtest
#include <string>
#include "gtest/gtest.h"
#include "qpragma.h"
#include "qpragma/hhl.hpp"

using namespace qpragma;


double ERROR_RATIO = 0.1;


/* HHL*/

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

TEST(HHL, m2x2_usecase2) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 I + 0.25 Z */
    hhl::observables::PauliTerm<SIZE> term1 { "I", 0.5 };
    hhl::observables::PauliTerm<SIZE> term2 { "Z", 0.25};
    std::vector term_vect { term1, term2 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0> + |1>
        std::array<double, 1UL << SIZE> init_array {1., 1.};
        qpragma::hhl::utils::normalize<SIZE>(init_array);

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0.10, 0.90};


        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m2x2_usecase3) {
    constexpr uint64_t SIZE = 1UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 I + 0.125 X */
    hhl::observables::PauliTerm<SIZE> term1 { "I", 0.5 };
    hhl::observables::PauliTerm<SIZE> term2 { "X", 0.125};
    std::vector term_vect { term1, term2 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    {
        // b = |0>
        std::array<double, 1UL << SIZE> init_array {1., 0.};

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0.94, 0.06};


        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}


/* 2 qubits - 4x4 matrices */


TEST(HHL, m4x4_II) {
    constexpr uint64_t SIZE = 2UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 II */
    hhl::observables::PauliTerm<SIZE> term1 { "II", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    for (uint64_t val = 0UL ; val < (1UL << SIZE) ; ++val) {
        // b = |val>
        std::array<double, 1UL << SIZE> init_array {0., 0., 0., 0.};
        init_array[val] = 1.;

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 0., 0., 0.};
        expected_result[val] = 1.;

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}

TEST(HHL, m4x4_XX) {
    constexpr uint64_t SIZE = 2UL;
    constexpr uint64_t SIZE_C = 4UL;
    constexpr uint64_t NB_SHOTS = 100UL;

    /* A = 0.5 II */
    hhl::observables::PauliTerm<SIZE> term1 { "XX", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

    for (uint64_t val = 0UL ; val < (1UL << SIZE) ; ++val) {
        // b = |val>
        std::array<double, 1UL << SIZE> init_array {0., 0., 0., 0.};
        init_array[val] = 1.;

        // Expected result
        std::array<double, 1UL << SIZE> expected_result {0., 0., 0., 0.};
        expected_result[(1UL << SIZE) - val - 1UL] = 1.;

        // Sampling on HHL algorithm
        auto res = hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, NB_SHOTS);

        for (uint64_t idx = 0UL ; idx < (1UL << SIZE) ; ++idx) {
            ASSERT_NEAR((double) res[idx] / NB_SHOTS, expected_result[idx], ERROR_RATIO);
        }
    }
}