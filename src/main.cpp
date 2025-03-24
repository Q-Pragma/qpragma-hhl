#include <iostream>

#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;
using namespace qpragma::hhl::stateprep;


// template <uint64_t SIZE>
// using basic_hhl = qpragma::hybrid_hhl<SIZE, kp_tree<SIZE>, trotterization<SIZE>>;

// DEFINE_HHL_IMPLEMENTATION(basic_hhl, kp_tree<SIZE>, trotterization<SIZE>);


void test_state_prep() {
    const uint64_t SIZE = 4UL;
    const uint64_t NB_SHOTS = 100;
    //init_array = qpragma::hhl::utils::normalize<SIZE>(init_array);
    std::vector<double> res(1 << SIZE);
    for (int i = 0 ; i < NB_SHOTS ; ++i) {
        qpragma::quint_t<SIZE> qreg;
        std::cout << "shot " << i << std::endl;
        std::array<double, 1 << SIZE> init_array {0.,0.,1.,1., 1.,0.,0.,0., 0., 0., 0., 0., 0., 0., 0., 0.8};
        init_array = qpragma::hhl::utils::normalize<SIZE>(init_array);
        (kp_tree<SIZE>(init_array))(qreg);
        auto idx = qpragma::measure_and_reset(qreg);
        res[idx] += 1.;
    }
    for (int i = 0 ; i < 1 << SIZE ; ++i) {
        std::cout << i << " : " << res[i] << std::endl;
    }
}

int main() {
    const uint64_t SIZE = 3UL;
    PauliTerm<SIZE> term { "XII", 5.2 };
    std::vector term_vect { term };
    std::array<double, 1 << SIZE> init_array {0.2, 0.3, 0.4, 0.5};

    Observable<SIZE> obs { term_vect };
    qpragma::array<SIZE> qreg;

    //basic_hhl<SIZE>(kp_tree<SIZE>(init_array), trotterization<SIZE>(obs, 0.3))(qreg);
    test_state_prep();

    qpragma::hhl::basic_hhl(qreg, init_array, obs);
}
