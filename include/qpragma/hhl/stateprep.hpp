/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/stateprep.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 *
 * @copyright   Licensed to the Apache Software Foundation (ASF) under one
 *              or more contributor license agreements. See the NOTICE file
 *              distributed with this work for additional information
 *              regarding copyright ownership. The ASF licenses this file
 *              to you under the Apache License, Version 2.0 (the
 *              "License"); you may not use this file except in compliance
 *              with the License.  You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 *              Unless required by applicable law or agreed to in writing,
 *              software distributed under the License is distributed on an
 *              "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 *              KIND, either express or implied.  See the License for the
 *              specific language governing permissions and limitations
 *              under the License.
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

#include <iostream>

#include "qpragma.h"
#include "qpragma/hhl/utils.hpp"


namespace qpragma::hhl::stateprep {
    const double _TOL = 1e-8;

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
        // 1-qubit system
        if constexpr (SIZE == 1UL) {
            double angle = 2. * qpragma::hhl::utils::sign(init_array[1]) 
                               * acos(qpragma::hhl::utils::sign(init_array[0]) * init_array[0]);
            (RY(angle))(qreg);
        }
        // Multi-qubit system
        else {
            std::vector<double> tree_vect = qpragma::hhl::stateprep::get_tree_coeff<SIZE>(init_array);
            // First rotation on the first qubit
            double angle = 2 * acos(sqrt(tree_vect[1]));
            if (std::abs(angle) > _TOL) {
                (RY(angle))(qreg[SIZE - 1UL]);
            }

            // Iter until the leaves (excluded) are reached --> work only with proba
            uint64_t start_val;
            double left, right;
            for (uint64_t idx = 1 ; idx < SIZE - 1 ; ++idx) {
                start_val = (1 << (idx + 1)) - 1;

                for (uint64_t ctrl_val = 0 ; ctrl_val < (1 << idx) ; ++ctrl_val) {
                    left = tree_vect[start_val + 2 * ctrl_val];
                    right = tree_vect[start_val + 2 * ctrl_val + 1];
                    if (left + right > _TOL) {
                        angle = 2 * acos(sqrt(left / (left + right)));
                        if (std::abs(angle) > _TOL) {
                            #pragma quantum ctrl (qpragma::as_uint<SIZE>(qreg(SIZE-idx, SIZE-1)) == ctrl_val)
                            (RY(angle))(qreg[SIZE - 1UL - idx]);
                        }
                    }
                }
            }

            // Last iteration : take into account signs from init_array
            start_val = (1 << SIZE) - 1;
            for (uint64_t ctrl_val = 0 ; ctrl_val < (1 << (SIZE-1)) ; ++ctrl_val) {
                double sign_left = qpragma::hhl::utils::sign(init_array[2 * ctrl_val]);
                double sign_right = qpragma::hhl::utils::sign(init_array[2 * ctrl_val + 1]);
                left = tree_vect[start_val + 2 * ctrl_val];
                right = tree_vect[start_val + 2 * ctrl_val + 1];
                if (left + right > _TOL) {
                    angle = 2 * sign_right * acos(sign_left * sqrt(left / (left + right)));
                    if (std::abs(angle) > _TOL) {
                        #pragma quantum ctrl (qpragma::as_uint<SIZE>(qreg(1, SIZE -1)) == ctrl_val)
                        (RY(angle))(qreg[0]);
                    }
                }
            }
        }
    }
}   // qpragma::hhl::stateprep


// Include implementation
#include "qpragma/hhl/stateprep.ipp"

#endif  /* QPRAGMA_HHL_STATEPREF_HPP */
