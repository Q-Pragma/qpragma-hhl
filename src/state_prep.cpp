#include <cmath>

#include "qpragma.h"
#include "utils.cpp"


/**
 * TREE APPROACH STATE PREP
 * Perform a state preparation on a n-qubit quantum state
 * according to a given vector of size 2^n. This is done
 * following the tree approach provided in KP16
 */

// Get the coefficients of the tree from the array
template <uint64_t SIZE>
std::vector<double> get_tree_coeff(std::array<double, (1 << SIZE)> init_vect) {
    uint64_t nb_terms = (1 << (SIZE + 1)) - 1;
    std::vector<double> tree_vect(nb_terms);
    // init the leaves
    std::transform(init_vect.begin(), init_vect.end(),
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

#pragma quantum routine (std::array<double, (1 << SIZE)> init_vect)
template <uint64_t SIZE>
void tree_state_prep(const quint_t<SIZE> & qreg) {
    if constexpr (SIZE == 1) {
         double angle = 2. * sign(init_vect[2]) * acos(sign(init_vect[1] * sqrt(init_vect[1])));
         (RY(angle))(qreg);
    }
    else {
        std::vector<double> tree_vect = get_tree_coeff<SIZE>(init_vect);
        // First rotation on the first qubit
        double angle = 2 * acos(sqrt(tree_vect[1]));
        (RY(angle))(qreg[0]);
        // Iter until the leaves (excluded) are reached --> work only with proba
        for (uint64_t idx = 1 ; idx < SIZE - 1 ; ++idx) {
            uint64_t start_val = (1 << (idx + 1)) - 1;
            for (uint64_t ctrl_val = 0 ; ctrl_val < (1 << idx) ; ++ctrl_val) {
                angle = 2 * acos(sqrt(tree_vect[start_val + 2 * ctrl_val]));
                #pragma quantum ctrl (qreg(0, idx-1) == ctrl_val)
                (RY(angle))(qreg[idx]);
            }
        }
        // Last iteration : take into account signs
        for (uint64_t ctrl_val = 0 ; ctrl_val << (1 < SIZE) ; ++ctrl_val) {
            angle = 2 * sign(init_vect[2 * ctrl_val + 1]) * acos(
                sign(init_vect[2 * ctrl_val]) * sqrt(tree_vect[ctrl_val + 2 * ctrl_val])
            );
            #pragma quantum ctrl (qreg(0, SIZE-2) == ctrl_val)
            (RY(angle))(qreg[SIZE - 1]);
        }
    }
}


// Main function to call state preparation
template <uint64_t SIZE>
void state_prep(const quint_t<SIZE> & qreg, std::array<double, (1 << SIZE)> init_vect) {
    (tree_state_prep<SIZE>(init_vect))(qreg);
}



