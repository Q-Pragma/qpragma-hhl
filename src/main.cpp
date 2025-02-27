#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;
using namespace qpragma::hhl::stateprep;


template <uint64_t SIZE>
using basic_hhl = qpragma::hybrid_hhl<SIZE, decltype(kp_tree<5UL>({})), decltype(trotterization<SIZE>(Observable<SIZE>(), 0.3))>;


int main() {
    PauliTerm<5UL> term { "XXXII", 5.2 };
    std::vector term_vect { term };
    std::array<double, 1 << 5> init_array {1.};

    Observable<5UL> obs { term_vect };
    qpragma::array<5UL> qreg;

    basic_hhl<5UL>(kp_tree<5UL>(init_array), trotterization<5UL>(obs, 0.3))(qreg);
}
