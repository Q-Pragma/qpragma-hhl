/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
 * @brief
 *
 */

#ifndef QPRAGMA_HHL_HPP
#define QPRAGMA_HHL_HPP

#include "qpragma.h"
// Include ...


namespace qpragma {
    #pragma quantum routine(STATE_PREP state_prep, HAM_SIM simu)
    template <uint64_t SIZE, typename STATE_PREP, typename HAM_SIM>
    void hybrid_hhl(const qpragma::array<SIZE> & qreg) {
        state_prep(qreg);
        simu(qreg);
    }
}

#endif  /* QPRAGMA_HLL_HPP */
