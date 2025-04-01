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
    if (val <= (1 << (nb_bits - 1))) {
        return (double) val / (double) (1 << (nb_bits - 1));
    }
    else {
        return (double) val / (double) ( 1 << (nb_bits - 1)) - 2.;
    }
}
