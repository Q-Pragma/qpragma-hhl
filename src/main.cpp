#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;
using namespace qpragma::hhl::stateprep;


template <uint64_t SIZE>
using basic_hhl = qpragma::hybrid_hhl<SIZE, decltype(kp_tree<SIZE>({})), decltype(trotterization<SIZE>(Observable<SIZE>(), 0.3))>;


int main() {
    const uint64_t SIZE = 3UL;
    PauliTerm<SIZE> term { "XII", 5.2 };
    std::vector term_vect { term };
    std::array<double, 1 << SIZE> init_array {0.2, 0.3, 0.4, 0.5};

    Observable<SIZE> obs { term_vect };
    qpragma::array<SIZE> qreg;

    basic_hhl<SIZE>(kp_tree<SIZE>(init_array), trotterization<SIZE>(obs, 0.3))(qreg);
}
