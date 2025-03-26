#include <algorithm>

/* Normalize an array of coefficients */
template <uint64_t SIZE>
std::array<double, (1 << SIZE)> qpragma::hhl::utils::normalize(const std::array<double, (1 << SIZE)> & coeffs) {
    
    double norm = 0.;
    for (uint64_t i = 0UL ; i < (1 << SIZE) ; ++i) {
        norm += coeffs[i] * coeffs[i];
    }
    norm = sqrt(norm);
    
    std::array<double, (1 << SIZE)> res(coeffs);
    std::transform(res.cbegin(), res.cend(),
                   res.begin(),
                   [norm](double e) {return e/norm; });

    return res;
}

/* Get the mean of each bit of all the eigenvalues */
template <uint64_t SIZE_C>
std::array<double, SIZE_C> qpragma::hhl::utils::get_means(const std::vector<uint64_t>& eigenvals) {

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