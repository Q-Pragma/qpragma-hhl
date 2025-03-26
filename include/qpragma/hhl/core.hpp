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
        for (int64_t idx = SIZE_C - 1; idx >= 0; --idx) {
            for (uint64_t loop_step = 0UL ; loop_step < (1UL << idx) ; ++loop_step) {
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
            std::cout << utils::bin_to_double(SIZE_C, i) << " : " << is_eigen[i] << std::endl;
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
    #pragma quantum routine (double c, std::array<double, SIZE_C> means)
    template <uint64_t SIZE_C>
    void reduced_AQE(const quint_t<SIZE_C> & creg, const qbool & anc) {
        for (uint64_t val_c = 1UL ; val_c < (1 << SIZE_C) ; ++val_c) {
            // Only on the compatible values with eigenvalues
            if (is_compatible<SIZE_C>(val_c, means, SIZE_C)) {
                // Angle of the rotation RY
                double theta = 2. * acos(sqrt(1 - c*c /(val_c*val_c)));
                // Rotation on the ancilla controlled by the eigenvalue
                #pragma quantum ctrl (creg == val_c)
                (RY(theta))(anc);
            }
        }
    }

    /* Reduced version of the QFT */
    #pragma quantum routine (std::array<double, SIZE_C> means)
    template <uint64_t SIZE_C>
    void reduced_qft(const quint_t<SIZE_C> & creg) {
        for (uint64_t target = 0UL ; target < SIZE_C ; ++target) {
            // On 0 nothing happens here
            if (means[target] == 0.)
                continue;

            // Apply an H gate where the mean is not fixed
            if (means[target] != 1.)
                H(creg[target]);

            // controlled phase part of the QFT
            for (uint64_t control = target + 1UL ; control < SIZE_C ; ++control) {
                double angle = M_PI / (1 << (control - target));
                if (means[control] == 1.) {
                    // No control needed because control qubit always 1
                    (PH(angle))(creg[target]);
                }
                else if (means[control] != 0.) {
                    #pragma quantum ctrl (creg[control])
                    (PH(angle))(creg[target]);
                }
            }
        }
    }

    /* Reduced QPE */
    #pragma quantum routine (const HAM_SIM & simu, std::array<double, SIZE_C> means)
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    void reduced_QPE(const quint_t<SIZE> & qreg, const quint_t<SIZE_C> & creg) {
        for (int64_t idx = SIZE_C - 1UL ; idx >= 0 ; --idx) {
            // if qubit i in 0 --> nothing happens
            if (means[idx] == 0.)
                continue;

            // if qubit i in 1 --> X
            if (means[idx] == 1.)
                X(creg[idx]);
            
            // if qubit mean not fixed --> H
            else
                H(creg[idx]);
            
            #pragma quantum ctrl (creg[idx])
            simu(qreg);
        }
        // Call the reduced QFT
        (reduced_qft<SIZE_C>(means))(creg);
    }

    /* Reduced version of HHL */
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM, typename STATE_PREP>
    void reduced_HHL(const HAM_SIM & simu,
                     const STATE_PREP & state_prep,
                     std::vector<uint64_t> eigenvals,
                     double c,
                     quint_t<SIZE> & qreg) {
        qbool anc = 0;
        // Get the means of the bits of the eigenvalues
        std::array<double, SIZE_C> means = utils::get_means<SIZE_C>(eigenvals);
        
        quint_t<SIZE_C> creg;
        // Post selection on ancilla in 1 state
        do {
            reset(qreg);
            state_prep(qreg);
            {
                #pragma quantum compute
                {
                    reduced_QPE<SIZE, SIZE_C, HAM_SIM>(simu, means)(qreg, creg);
                }
        
                (reduced_AQE<SIZE_C>(c, means))(creg, anc);
            }
            // Automatically uncompute reduced_QPE
            
            reset(creg);
        } while (not measure_and_reset(anc));
    }


    template <uint64_t SIZE, uint64_t SIZE_C, typename STATE_PREP, typename HAM_SIM>
    void hybrid_hhl(
        qpragma::quint_t<SIZE> & qreg,
        const std::array<double, (1UL << SIZE)> & init,
        const qpragma::hhl::observables::Observable<SIZE> & observable
    ) {
        // Initialize the state preparation and hamiltonian simulation
        STATE_PREP state_prep { init };
        HAM_SIM simu { observable };

        double c = 3.;  // TODO: trouver la valeur de c

        // Hybrid quantum-classical sampling of eigenvalues
        std::vector<uint64_t> eigenvals = get_eigenvals<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep);
        
        std::cout << "Start printing eigenvals :" << std::endl;
        for (uint64_t i = 0UL; i < eigenvals.size() ; ++i) {
            std::cout << i << " : " << utils::bin_to_double(SIZE_C, eigenvals[i]) << std::endl;
        }
        std::cout << "End printing eigenvals" << std::endl;
        
        // Call to reduced HHL
        //reduced_HHL<SIZE, SIZE_C, HAM_SIM, STATE_PREP>(simu, state_prep, eigenvals, c, qreg);
    }
}


#endif  /* QAT_ */

