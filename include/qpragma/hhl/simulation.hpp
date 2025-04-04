/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/simulation.hpp
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
 * Implementation of Hamiltonian simulation routines, using Q-Pragma
 */

#ifndef QPRAGMA_HHL_SIMULATION_HPP
#define QPRAGMA_HHL_SIMULATION_HPP

#include <ranges>
#include <cmath>

#include "qpragma.h"
#include "qpragma/hhl/observables.hpp"


namespace qpragma::hhl::simulation {
    /**
     * Given a Pauli term P, apply exp(-i P)
     * to a quantum register
     */
    #pragma quantum routine(observables::PauliTerm<SIZE> term)
    template <uint64_t SIZE>
    void apply_term(const qpragma::array<SIZE> & qreg) {
        using namespace qpragma;
        using qpragma::hhl::observables::pauli_op;

        // Apply a global phase if PauliTerm is (coeff, "III...I")
        // This global phase is useful when this routine is controlled
        if (term.is_identity()) {
            qbool ancilla = true;
            PH(- M_PI * term.coeff())(ancilla);
        }

        // Apply term
        else {
            // Get index of the first qubit having an operator no equal to I
            ssize_t first_qubit = -1;

            // Apply base change
            #pragma quantum compute
            {
                for (int64_t idx = SIZE - 1; idx >= 0; --idx) {
                    switch (term[idx]) {
                    case pauli_op::I:
                        continue;
                    case pauli_op::X:
                        H(qreg[idx]);
                        break;
                    case pauli_op::Y:
                        RX(-M_PI / 2.)(qreg[idx]);
                        break;
                    default:
                        break;
                    };
                    // Update "first_qbit" or apply CNOT
                    if (first_qubit == -1)
                        first_qubit = static_cast<ssize_t>(idx);
                    else
                        CNOT(qreg[idx], qreg[first_qubit]);
                }
            }

            RZ(2. * M_PI * term.coeff())(qreg[first_qubit]);
        }
    }


    /**
     * Perform hamiltonian simulation using a trotterization
     */
    #pragma quantum routine(observables::Observable<SIZE> observable, double eps = 0.1)
    template <uint64_t SIZE>
    void trotterization(const qpragma::array<SIZE> & qreg) {
        uint64_t n_trotter = 1. / sqrt(eps);

        observable *= -1. / static_cast<double>(n_trotter);

        for (uint64_t t_step = 0UL; t_step < n_trotter; ++t_step) {
            for (auto term: observable) {
                (apply_term<SIZE>(term))(qreg);
            }
        }
    }
}   // qpragma::hhl::simulation

#endif  /* QPRAGMA_HHL_SIMULATION_HPP */
