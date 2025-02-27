#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;


template <uint64_t SIZE>
using basic_hhl = qpragma::hybrid_hhl<SIZE, decltype(tree_state_prep<5UL>({})), decltype(trotterization<SIZE>(Observable<SIZE>(), 0.3))>;


int main() {
    PauliTerm<5UL> term { "XXXII", 5.2 };
    std::vector term_vect { term };
    std::array<double, 1 << 5> init_array {1.};

    Observable<5UL> obs { term_vect };
    qpragma::array<5UL> qreg;

    basic_hhl<5UL>(tree_state_prep<5UL>(init_array), trotterization<5UL>(obs, 0.3))(qreg);
}
