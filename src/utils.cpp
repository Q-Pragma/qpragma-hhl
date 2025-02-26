#include <cstdint>
#include <array>
#include <cmath>

/* Get the sign of a double*/
double sign(const double & d) {
    return d > 0. ? 1. : -1.;
}

/* Normalize an array of coefficients */
template <uint64_t SIZE>
std::array<double, (1 << SIZE)> normalize(std::array<double, (1 << SIZE)> coeffs) {
    double norm = 0.;
    std::array<double, (1 << SIZE)> res(coeffs);
    for (int i = 0 ; i < coeffs.size() ; ++i) {
        norm += coeffs[i] * coeffs[i];
    }
    norm = sqrt(norm);
    std::transform(res.cbegin(), res.cend(),
                   res.begin(),
                   [norm](double e) {return e/norm; });

    return res;
}

/* Convert the binary value to a double flotting point value */
double bin_to_double(uint64_t nb_bits, uint64_t val) {
    if (val == 0) {
        return 1.;  // 0.000 and 1.000 equivalent but only 1. can be eigenvalue
    }
    else {
        return (double) val / (double) (1 << nb_bits);
    }
}

/* Get the mean of each bit of all the eigenvalues */
template <uint64_t SIZEC>
std::array<double, SIZEC> get_means(std::vector<uint64_t> eigenvals) {
    std::array<double, SIZEC> means;
    for (int i = 0 ; i < SIZEC ; ++i) {
        means[i] = 0.;
        // Compute the mean on bit i
        for (uint64_t val : eigenvals) {
            means[i] += (double) ((val >> i) & 1);
        }
        means[i] /= (double) eigenvals.size();
    }
    return means;
}
