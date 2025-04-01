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
        init_array = normalize<SIZE>(init_array);
        (kp_tree<SIZE>(init_array))(qreg);
        bool idx = qpragma::measure_and_reset(qreg);
        res[idx] += 1.;
    }
}


int main() {
    const uint64_t SIZE = 2UL;
    const uint64_t SIZE_C = 4UL;
    PauliTerm<SIZE> term1 { "II", 0. };
    PauliTerm<SIZE> term2 { "XI", 0.5 };
    PauliTerm<SIZE> term3 { "IX", 0.25 };
    std::vector term_vect { term1, term2, term3 };
    std::array<double, 1 << SIZE> init_array {1., 0., 0., 0};

    Observable<SIZE> obs { term_vect };
    qpragma::quint_t<SIZE> qreg;

    qpragma::hhl::basic_hhl<SIZE, SIZE_C>(qreg, init_array, obs);
}

