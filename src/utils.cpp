#include <cstdint>
#include <array>
#include <cmath>

#include "qpragma/hhl/utils.hpp"

/* Get the sign of a double*/
double qpragma::hhl::utils::sign(const double & d) {
    return d > 0. ? 1. : -1.;
}

/* Convert the binary value to a double flotting point value */
double qpragma::hhl::utils::bin_to_double(uint64_t nb_bits, uint64_t val) {
    if (val == 0) {
        return 1.;  // 0.000 and 1.000 equivalent but only 1. can be eigenvalue
    }
    else {
        return (double) val / (double) (1 << nb_bits);
    }
}
