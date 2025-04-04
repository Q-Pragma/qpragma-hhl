/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl.hpp
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
 *
 */

#ifndef QPRAGMA_HHL_HPP
#define QPRAGMA_HHL_HPP

#include "qpragma.h"
#include "qpragma/hhl/core.hpp"
#include "qpragma/hhl/observables.hpp"
#include "qpragma/hhl/simulation.hpp"
#include "qpragma/hhl/stateprep.hpp"

#define DEFINE_HHL_IMPLEMENTATION(name, state_prep_t, simu_t) \
    template <uint64_t SIZE, uint64_t SIZE_C> \
    void name ( \
        qpragma::quint_t<SIZE> & qreg, \
        const std::array<double, (1UL << SIZE)> & init, \
        const qpragma::hhl::observables::Observable<SIZE> & observable, \
        double eps = 0.1 \
    ) { \
        qpragma::hhl::hybrid_hhl<SIZE, SIZE_C, decltype(state_prep_t{ init }), decltype(simu_t{ observable })>( \
            qreg, init, observable, eps \
        ); \
    }

#define DEFINE_HHL_IMPLEMENTATION_WITH_SAMPLING(name, state_prep_t, simu_t) \
    template <uint64_t SIZE, uint64_t SIZE_C> \
    std::array<uint64_t, 1UL << SIZE> name ( \
        const std::array<double, (1UL << SIZE)> & init, \
        const qpragma::hhl::observables::Observable<SIZE> & observable, \
        uint64_t nb_shots = 1UL, \
        double eps = 0.1 \
    ) { \
        return qpragma::hhl::hybrid_hhl_with_sampling<SIZE, SIZE_C, decltype(state_prep_t{ init }), decltype(simu_t{ observable })>( \
            init, observable, nb_shots, eps \
        ); \
    }


namespace qpragma::hhl {
    DEFINE_HHL_IMPLEMENTATION(basic_hhl, stateprep::kp_tree<SIZE>, simulation::trotterization<SIZE>)
    DEFINE_HHL_IMPLEMENTATION_WITH_SAMPLING(test_hhl, stateprep::kp_tree<SIZE>, simulation::trotterization<SIZE>)
}

#endif  /* QPRAGMA_HLL_HPP */
