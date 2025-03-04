/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/observables.hpp
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
