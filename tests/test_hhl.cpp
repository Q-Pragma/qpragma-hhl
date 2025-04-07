// Unitary test using gtest
#include <string>
#include "gtest/gtest.h"
#include "qpragma.h"
#include "qpragma/hhl.hpp"

using namespace qpragma;


double ERROR_RATIO = 0.1;


TEST(HHL, m1x1_I) {
    const uint64_t SIZE = 1UL;
    const uint64_t SIZE_C = 4UL;
    const uint64_t NB_SHOTS = 100UL;

    /* A = I */
    hhl::observables::PauliTerm<SIZE> term1 { "I", 0.5 };
    std::vector term_vect { term1 };
    hhl::observables::Observable<SIZE> obs { term_vect };

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