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


namespace qpragma::hhl {
    /**
     * Controlled hamiltonian simulation
     * This function performs a QPE (without QFT at the end)
     */
    #pragma quantum routine (const HAM_SIM & simu)
    template <uint64_t SIZE, uint64_t SIZE_C, typename HAM_SIM>
    void controlled_simu(const array<SIZE> & reg, const array<SIZE_C> & creg) {
        for (uint64_t idx = 0UL; idx < SIZE_C; ++idx) {
            #pragma quantum ctrl(creg[SIZE_C - idx - 1UL])
            simu(reg);  // TODO: For loop
        }
    }


    template <uint64_t SIZE, typename STATE_PREP, typename HAM_SIM>
    void hybrid_hhl(
        const qpragma::array<SIZE> & qreg,
        const std::array<double, (1UL << SIZE)> & init,
        const qpragma::hhl::observables::Observable<SIZE> & observable
    ) {
        STATE_PREP state_prep { init };
        HAM_SIM simu { observable };

        state_prep(qreg);
        // simu(qreg);
        //

        array<SIZE> creg;
        (controlled_simu<SIZE, SIZE, HAM_SIM>(simu))(qreg, creg);
    }
}


#endif  /* QAT_ */

