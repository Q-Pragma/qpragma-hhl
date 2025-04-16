/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/core.hpp
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
 */

#ifndef QPRAGMA_HHL_CORE_HPP
#define QPRAGMA_HHL_CORE_HPP

#include "qpragma.h"
#include "qpragma/hhl/observables.hpp"
#include "qpragma/hhl/utils.hpp"


namespace qpragma::hhl {
    const uint64_t _NB_SHOTS = 100UL;

    /**
     * Controlled hamiltonian simulation
     * This function performs a QPE (without QFT at the end)
     */
    #pragma quantum routine (const HAM_SIM & simu)
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    void controlled_simu(const quint_t<SIZE> & qreg, const quint_t<SIZE_C> & creg) {
        constexpr int64_t signed_size = static_cast<int64_t>(SIZE_C);
        for (int64_t idx = signed_size - 1L ; idx >= 0L ; --idx) {
            for (int64_t loop_step = 0L ; loop_step < (1L << idx) ; ++loop_step) {
                #pragma quantum ctrl(creg[idx])
                simu(qreg);
            }
        }
    }

    /* Quantum Phase Estimation algorithm */
    #pragma quantum routine (const HAM_SIM & simu)
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    void QPE(const quint_t<SIZE> & qreg, const quint_t<SIZE_C> & creg) {
        wall::H<SIZE_C>(creg);
        (controlled_simu<SIZE, SIZE_C, HAM_SIM>(simu))(qreg, creg);
        qft<SIZE_C>.dag(creg);
    }

    /* Quantum Phase Estimation with a measurement of eigenvalue */
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    uint64_t QPEA(const HAM_SIM & simu, const quint_t<SIZE> & qreg) {
        quint_t<SIZE_C> creg = 0UL;
        (QPE<SIZE, SIZE_C, HAM_SIM>(simu))(qreg, creg);
        return measure_and_reset(creg);
    }

    /* Get an estimation of all the eigenvalues by sampling on QPE */
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM, typename STATE_PREP>
    std::vector<uint64_t> get_eigenvals(const HAM_SIM & simu, const STATE_PREP & state_prep) {
        // Initialize an array to check if a value is an eigenvalue
        std::array<uint64_t, 1 << SIZE_C> is_eigen;
        is_eigen.fill(0UL);

        quint_t<SIZE> qreg;
        // Sample on QPE
        for (uint64_t i = 0UL; i < _NB_SHOTS ; ++i) {
            state_prep(qreg);
            uint64_t res = QPEA<SIZE, SIZE_C, HAM_SIM>(simu, qreg);
            ++is_eigen[res];
            reset(qreg);
        }

        // Create the vector containing the eigenvalues
        std::vector<uint64_t> eigenvals;
        for (uint64_t i = 0UL; i < (1 << SIZE_C) ; ++i) {
            if (is_eigen[i])
                eigenvals.push_back(i);
        }

        return eigenvals;
    }

    /* Check if a value lambda is compatible with the means observed on eigenvalues */
    template <uint64_t SIZE_C>
    bool is_compatible(uint64_t lambda, std::array<double, SIZE_C> means, const uint64_t SIZE) {
        for (int i = 0 ; i < SIZE ; ++i) {
            if (means[i] == 0. or means[i] == 1.) {
                if (((lambda >> i) & 1) != means[i])
                    return false;  // Is not compatible with the eigenvalues
            }
        }
        return true;  // Is compatible
    }

    /* Reduced version of the AQE */
    #pragma quantum routine (double c, std::vector<uint64_t> eigenvals)
    template <uint64_t SIZE_C>
    void reduced_AQE(const quint_t<SIZE_C> & creg, const qbool & anc) {
        for (auto val_c : eigenvals) {
            // Angle of the rotation RY
            double val_c_d = utils::bin_to_double(SIZE_C, val_c);
            double theta = 2 * utils::sign(val_c_d) * acos(sqrt(1 - c*c / (val_c_d * val_c_d)));
            // Rotation on the ancilla controlled by the eigenvalue
            #pragma quantum ctrl (creg == val_c)
            (RY(theta))(anc);
        }
    }

    /* Reduced version of the QFT */
    #pragma quantum routine (std::array<double, SIZE_C> means)
    template <uint64_t SIZE_C>
    void reduced_qft(const quint_t<SIZE_C> & creg) {
        constexpr int64_t signed_size = static_cast<int64_t>(SIZE_C);
        for (int64_t target = signed_size - 1L ; target >= 0L ; --target) {
            // On 0 nothing happens here
            if (means[target] == 0.)
                continue;

            // Apply an H gate where the mean is not fixed
            if (means[target] != 1.)
                H(creg[target]);

            // controlled phase part of the QFT
            for (int64_t control = target - 1L ; control >= 0L ; --control) {
                double angle = M_PI / (1 << (target - control));
                if (means[control] == 1.) {
                    // No control needed because control qubit always 1
                    (PH(angle))(creg[target]);
                }
                else if (means[control] != 0.) {
                    (PH(angle)).ctrl(creg[control], creg[target]);
                }
            }
        }
        reverse_order<SIZE_C>(creg);
    }

    /* Reduced QPE */
    #pragma quantum routine (const HAM_SIM & simu, std::array<double, SIZE_C> means)
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    void reduced_QPE(const quint_t<SIZE> & qreg, const quint_t<SIZE_C> & creg) {
        constexpr int64_t signed_size = static_cast<int64_t>(SIZE_C);
        for (int64_t idx = signed_size - 1L ; idx >= 0L ; --idx) {
            // if qubit i in 0 --> nothing happens
            if (means[idx] == 0.)
                continue;

            // if qubit i in 1 --> X
            if (means[idx] == 1.)
                X(creg[idx]);
            
            // if qubit mean not fixed --> H
            else
                H(creg[idx]);
            
            for (int64_t loop_step = 0L ; loop_step < (1L << idx) ; ++loop_step) {
                #pragma quantum ctrl (creg[idx])
                simu(qreg);
            }
        }
        // Call the reduced QFT
        (reduced_qft<SIZE_C>(means)).dag(creg);
    }

    /* Reduced version of HHL */
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM, typename STATE_PREP>
    void reduced_HHL(const HAM_SIM & simu,
                     const STATE_PREP & state_prep,
                     std::vector<uint64_t> eigenvals,
                     double c,
                     quint_t<SIZE> & qreg) {
        
        qbool anc = 0;
        // Post selection on ancilla in 1 state
        do {
            reset(qreg);
            state_prep(qreg);
            {
                quint_t<SIZE_C> creg = 0UL;

                #pragma quantum compute
                {
                    (QPE<SIZE, SIZE_C, HAM_SIM>(simu))(qreg, creg);
                }
        
                (reduced_AQE<SIZE_C>(c, eigenvals))(creg, anc);

                // Automatically uncompute QPE
            }
        } while (not measure_and_reset(anc));
    }


    template <uint64_t SIZE, uint64_t SIZE_C, typename STATE_PREP, typename HAM_SIM>
    void hybrid_hhl(
        qpragma::quint_t<SIZE> & qreg,
        const std::array<double, (1UL << SIZE)> & init,
        const qpragma::hhl::observables::Observable<SIZE> & observable,
        double eps = 0.1
    ) {
        // Initialize the state preparation and hamiltonian simulation
        STATE_PREP state_prep { init };
        HAM_SIM simu { observable, eps };

        // Hybrid quantum-classical sampling of eigenvalues
        std::vector<uint64_t> eigenvals = get_eigenvals<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep);
        
        // Find the smallest eigenvalue (in abs)
        double c = 1.;
        for (uint64_t i = 0UL; i < eigenvals.size() ; ++i) {
            double abs_eig = fabs(utils::bin_to_double(SIZE_C, eigenvals[i]));
            if (abs_eig < c)
                c = abs_eig;
        }

        reduced_HHL<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep, eigenvals, c, qreg);
    }

    template <uint64_t SIZE, uint64_t SIZE_C, typename STATE_PREP, typename HAM_SIM>
    std::array<uint64_t, 1UL << SIZE> hybrid_hhl_with_sampling(
        const std::array<double, (1UL << SIZE)> & init,
        const qpragma::hhl::observables::Observable<SIZE> & observable,
        uint64_t nb_shots = 1UL,
        double eps = 0.1
    ) {
        // Initialize the state preparation and hamiltonian simulation
        STATE_PREP state_prep { init };
        HAM_SIM simu { observable, eps };

        // Hybrid quantum-classical sampling of eigenvalues
        std::vector<uint64_t> eigenvals = get_eigenvals<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep);

        // Find the smallest eigenvalue (in abs)
        double c = 1.;
        for (uint64_t i = 0UL; i < eigenvals.size() ; ++i) {
            double abs_eig = fabs(utils::bin_to_double(SIZE_C, eigenvals[i]));
            if (abs_eig < c)
                c = abs_eig;
        }

        quint_t<SIZE> qreg;
        std::array<uint64_t, 1UL << SIZE> res;
        res.fill(0UL);
        for (uint64_t step = 0UL; step < nb_shots ; ++step) {
            reduced_HHL<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep, eigenvals, c, qreg);
            ++res[measure_and_reset(qreg)];
        }
        return res;
    }
}

#endif  /* QAT_ */

