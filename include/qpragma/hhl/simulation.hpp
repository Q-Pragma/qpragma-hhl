/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/simulation.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
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
            PH(term.coeff())(ancilla);
        }

        // Apply term
        else {
            std::vector<uint64_t> cnot_keep;  // idx where pauli_op != I

            // Apply base change
            #pragma quantum compute
            {
                for (uint64_t idx = 0UL; idx < SIZE; ++idx) {
                    switch (term[idx]) {
                    case pauli_op::I:
                        continue;
                    case pauli_op::X:
                        H(qreg[idx]);
                        cnot_keep.push_back(idx);
                        break;
                    case pauli_op::Y:
                        RX(-M_PI / 2.)(qreg[idx]);
                        cnot_keep.push_back(idx);
                        break;
                    default:
                        cnot_keep.push_back(idx);
                        break;
                    };
                }
            }

            #pragma quantum compute
            {
                for (uint64_t idx = 0 ; idx < cnot_keep.size() - 1UL ; ++idx) {
                    CNOT(qreg[cnot_keep[idx]], qreg[cnot_keep[idx + 1]]);
                }
            }
            RZ(term.coeff())(qreg[cnot_keep.back()]);
        }
    }


    /**
     * Perform hamiltonian simulation using a trotterization
     */
    #pragma quantum routine(observables::Observable<SIZE> observable, double esp)
    template <uint64_t SIZE>
    void trotterization(const qpragma::array<SIZE> & qreg) {
        uint64_t n_trotter = 10UL;

        observable *= 1. / static_cast<double>(n_trotter);

        for (uint64_t t_step = 0UL; t_step < n_trotter; ++t_step) {
            for (auto term: observable) {
                (apply_term<SIZE>(term))(qreg);
            }
        }
    }
}   // qpragma::hhl::simulation

#endif  /* QPRAGMA_HHL_SIMULATION_HPP */
