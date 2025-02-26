/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/observables.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
 *
 * @brief
 * Implementation of observables using Q-Pragma
 */

#ifndef QPRAGMA_HHL_OBSERVABLES_HPP
#define QPRAGMA_HHL_OBSERVABLES_HPP

#include <span>
#include <array>
#include <vector>
#include <cstdint>


namespace qpragma::hhl::observables {
    /**
     * TODO
     */
    enum class pauli_op { I='I', X='X', Y='Y', Z='Z' };


    /**
     * TODO
     */ 
    template <uint64_t SIZE>
    class PauliTerm {
    public:
        // Pauli operator

    private:
        std::array<pauli_op, SIZE> m_pauli_string;
        double m_coeff;

    public:
        PauliTerm(std::span<const pauli_op, SIZE> /* pauli_string */, double /* coeff */);
    };


    /**
     * TODO
     */
    template <uint64_t SIZE>
    class Observable {
    private:
        std::vector<PauliTerm<SIZE>> m_pauli_terms;

    public:
        // Constructor
        Observable(const std::vector<PauliTerm<SIZE>> & /* pauli_terms */);

        // Getters
        // std::vector<double> get_coeffs() const;
        // std::vector<std::array<pauli_op, SIZE>> get_strings() const;
    };
}   // qpragma::hhl::observables

#endif  /* QPRAGMA_HHL_OBSERVABLES_HPP */
