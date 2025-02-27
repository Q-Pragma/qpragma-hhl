/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/stateprep.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
 *
 * @brief
 * Implementation of state-preparation routines, using Q-Pragma
 */

#ifndef QPRAGMA_HHL_STATEPREP_HPP
#define QPRAGMA_HHL_STATEPREP_HPP


namespace qpragma::hhl::stateprep {

    /**
     * TREE APPROACH STATE PREP
     * Perform a state preparation on a n-qubit quantum state
     * according to a given vector of size 2^n. This is done
     * following the tree approach provided in KP16
     */

    // Get the tree coefficients from a given array of size 2^SIZE
    template <uint64_t SIZE>
    std::vector<double> get_tree_coeff(std::array<double, (1 << SIZE)> /*init_array*/);

}   // qpragma::hhl::stateprep

#include "stateprep.ipp"

#endif  /* QPRAGMA_HHL_STATEPREF_HPP */
