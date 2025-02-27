// Tree based state preparation

template <uint64_t SIZE>
inline std::vector<double> qpragma::hhl::stateprep::get_tree_coeff(std::array<double, (1 << SIZE)> init_array) {
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

#pragma quantum routine (std::array<double, (1 << SIZE)> init_array)
template <uint64_t SIZE>
void tree_state_prep(const quint_t<SIZE> & qreg) {
    if constexpr (SIZE == 1) {
         double angle = 2. * sign(init_array[2]) * acos(sign(init_array[1] * sqrt(init_array[1])));
         (RY(angle))(qreg);
    }
    else {
        std::vector<double> tree_vect = qpragma::hhl::stateprep::get_tree_coeff<SIZE>(init_array);
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
            angle = 2 * sign(init_array[2 * ctrl_val + 1]) * acos(
                sign(init_array[2 * ctrl_val]) * sqrt(tree_vect[ctrl_val + 2 * ctrl_val])
            );
            #pragma quantum ctrl (qreg(0, SIZE-2) == ctrl_val)
            (RY(angle))(qreg[SIZE - 1]);
        }
    }
}



