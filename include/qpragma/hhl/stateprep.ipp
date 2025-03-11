#include <cstdint>
#include <array>
#include <cmath>

#include "qpragma/hhl/utils.hpp"

// Tree based state preparation

template <uint64_t SIZE>
inline std::vector<double> qpragma::hhl::stateprep::get_tree_coeff(const std::array<double, (1UL << SIZE)> & init_array) {
    uint64_t nb_terms = (1 << (SIZE + 1)) - 1;
    std::vector<double> tree_vect(nb_terms);
    // init the leaves
    std::transform(init_array.begin(), init_array.end(),
                   tree_vect.begin() + (1 << SIZE) - 1,
                   [](double val){ return val * val; });

    // iter until the root is reached
    for (int64_t l = SIZE - 1 ; l >= 0 ; --l) {
        uint64_t start_copy = (1 << (l+1)) - 1;
        uint64_t start_store = (1 << l) - 1;
        for (uint64_t idx = 0 ; idx < (1 << l) ; ++idx) {
            tree_vect[start_store + idx] = tree_vect[start_copy + 2*idx] + tree_vect[start_copy + 2*idx + 1];
        }
    }

    return tree_vect;
}
