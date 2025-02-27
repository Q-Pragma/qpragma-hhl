#include <algorithm>

/* Normalize an array of coefficients */
template <uint64_t SIZE>
std::array<double, (1 << SIZE)> qpragma::hhl::utils::normalize(std::array<double, (1 << SIZE)> coeffs) {
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

/* Get the mean of each bit of all the eigenvalues */
template <uint64_t SIZEC>
std::array<double, SIZEC> qpragma::hhl::utils::get_means(std::vector<uint64_t> eigenvals) {
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