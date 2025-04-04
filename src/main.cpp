

#include <iostream>

#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;
using namespace qpragma::hhl::stateprep;
using namespace qpragma::hhl::utils;


// template <uint64_t SIZE>
// using basic_hhl = qpragma::hybrid_hhl<SIZE, kp_tree<SIZE>, trotterization<SIZE>>;

// DEFINE_HHL_IMPLEMENTATION(basic_hhl, kp_tree<SIZE>, trotterization<SIZE>);

template<uint64_t SIZE>
void print_state(const std::array<double, 1 << SIZE> & arr) {
    for (uint64_t i = 0UL; i < 1 << SIZE ; ++i) {
        std::cout << i << " : " << arr[i] << std::endl;
    }
}

template<uint64_t SIZE, uint64_t NB_SHOTS>
void test_state_prep(std::array<double, 1 << SIZE> & init_array) {

    std::vector<double> res(1 << SIZE);
    for (uint64_t i = 0UL; i < NB_SHOTS ; ++i) {
        qpragma::quint_t<SIZE> qreg;
        normalize<SIZE>(init_array);
        (kp_tree<SIZE>(init_array))(qreg);
        bool idx = qpragma::measure_and_reset(qreg);
        res[idx] += 1.;
    }
}


int main() {
    const uint64_t SIZE = 2UL;
    const uint64_t SIZE_C = 4UL;
    PauliTerm<SIZE> term1 { "IX", 0.5 };
    PauliTerm<SIZE> term2 {"XI" , 0.25};
    std::vector term_vect { term1, term2 };
    std::array<double, 1 << SIZE> init_array {1., 0., 0., 0.};

    Observable<SIZE> obs { term_vect };

    auto res = qpragma::hhl::test_hhl<SIZE, SIZE_C>(init_array, obs, 100UL);
    for (uint64_t idx =0 ; idx < (1 << SIZE) ; ++idx) {
        std::cout << idx << " " << res[idx] << std::endl;
    }
}

