#include <algorithm>
#include <cmath>


/* Normalize an array of coefficients */
template <uint64_t SIZE>
void qpragma::hhl::utils::normalize(std::array<double, (1 << SIZE)> & coeffs) {
    // Compute norm
    double norm = 0.;

    for (auto coeff: coeffs) {
        norm += coeff * coeff;
    }

    norm = std::sqrt(norm);

    // Normalize vector
    std::for_each(coeffs.begin(), coeffs.end(), [norm](double & amp) {return amp /= norm; });
}


/* Get the mean of each bit of all the eigenvalues */
template <uint64_t SIZE_C>
std::array<double, SIZE_C> qpragma::hhl::utils::get_means(const std::vector<uint64_t> & eigenvals) {
    size_t nb_eigenvals = eigenvals.size();
    std::array<double, SIZE_C> means{};

    for (uint64_t i = 0UL ; i < SIZE_C ; ++i) {
        double mean = 0.;

        // Compute the mean on bit i
        for (uint64_t val : eigenvals) {
            mean += (val >> i) & 1;
        }

        means[i] = mean / nb_eigenvals;
    }

    return means;
}
