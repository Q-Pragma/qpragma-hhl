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
#include <algorithm>


namespace qpragma::hhl::observables {
    /**
     * Enumeration describing a Pauli operator
     * A Pauli operator can be easily defined from a char
     */
    enum class pauli_op { I='I', X='X', Y='Y', Z='Z' };


    /**
     * Structure describing a Pauli String
     * A Pauli string is a sequence (of fixed size) of Pauli operators
     */ 
    template <uint64_t SIZE>
    class PauliTerm {
    private:
        std::array<pauli_op, SIZE> m_pauli_string;
        double m_coeff;

    public:
        PauliTerm(std::span<const char> /* pauli_string */, double /* coeff */);
        PauliTerm(const std::array<pauli_op, SIZE> & /* pauli_string */, double /* coeff */);

        // Methods
        double coeff() const;
        std::array<pauli_op, SIZE> pauli_string() const;

        std::array<pauli_op, SIZE>::const_iterator begin() const;
        std::array<pauli_op, SIZE>::const_iterator end() const;

        bool is_identity() const;

        pauli_op operator[](uint64_t /* idx */) const;
    };


    /**
     * Structure describing an observable
     * An observable is a set of Pauli strings
     */
    template <uint64_t SIZE>
    class Observable {
    private:
        std::vector<PauliTerm<SIZE>> m_pauli_terms;

    public:
        // Constructor
        Observable() = default;
        Observable(const std::vector<PauliTerm<SIZE>> & /* pauli_terms */);

        // Get number of terms
        uint64_t size() const;

        std::vector<PauliTerm<SIZE>>::const_iterator begin() const;
        std::vector<PauliTerm<SIZE>>::const_iterator end() const;

        Observable<SIZE>& operator*=(double /* constant */);
    };
}   // qpragma::hhl::observables


// Operator (double * Observable<SIZE>)
template <uint64_t SIZE>
qpragma::hhl::observables::Observable<SIZE> operator*(
    double /* constant */, const qpragma::hhl::observables::Observable<SIZE> & /* obs */
);


#include "qpragma/hhl/observables.ipp"

#endif  /* QPRAGMA_HHL_OBSERVABLES_HPP */
