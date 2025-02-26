#include "qpragma/hhl.hpp"
#include "qpragma/hhl/simulation.hpp"

using namespace qpragma::hhl::observables;
using namespace qpragma::hhl::simulation;


template <uint64_t SIZE>
using basic_hhl = qpragma::hybrid_hhl<SIZE, decltype(trotterization<SIZE>(Observable<SIZE>(), 0.3)), decltype(trotterization<SIZE>(Observable<SIZE>(), 0.3))>;


int main() {
    PauliTerm<5UL> term { "XXXII", 5.2 };
    std::vector term_vect { term };

    Observable<5UL> obs { term_vect };
    qpragma::array<5UL> qreg;

    basic_hhl<5UL>(trotterization<5UL>(obs, 0.3), trotterization<5UL>(obs, 0.3))(qreg);
}
