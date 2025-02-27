/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/stateprep.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
 *
 * @brief
 * Implementation of state-preparation routines, using Q-Pragma
 */

#ifndef QPRAGMA_HHL_STATEPREP_HPP
#define QPRAGMA_HHL_STATEPREP_HPP

#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

#include "qpragma.h"


namespace qpragma::hhl::stateprep {

    /**
     * TREE APPROACH STATE PREP
     * Perform a state preparation on a n-qubit quantum state
     * according to a given vector of size 2^n. This is done
     * following the tree approach provided in KP16
     */

    // Get the tree coefficients from a given array of size 2^SIZE
    template <uint64_t SIZE>
    inline std::vector<double> get_tree_coeff(const std::array<double, (1UL << SIZE)> & /* init_array */);


    #pragma quantum routine (std::array<double, (1 << SIZE)> init_array)
    template <uint64_t SIZE>
    void kp_tree(const array<SIZE> & qreg) {
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
}   // qpragma::hhl::stateprep


// Include implementation
#include "qpragma/hhl/stateprep.ipp"

#endif  /* QPRAGMA_HHL_STATEPREF_HLL */
