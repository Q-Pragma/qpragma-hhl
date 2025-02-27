/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/utils.hpp
 * @authors     Arnaud GAZDA <arnaud.gazda@eviden.com>
 * @copyright   2025 Bull S.A.S. - All rights reserved.
 *              This is not Free or Open Source software.
 *              Please contact BULL SAS for details about its license.
 *              Bull - Rue Jean Jaurès - B.P. 68 - 78340 Les Clayes-sous-Bois
 *
 * @brief
 * Util routines for Q-Pragma HHL
 */

 #ifndef QPRAGMA_HHL_UTILS_HPP
 #define QPRAGMA_HHL_UTILS_HPP
 
 #include <array>
 #include <vector>
 #include <cstdint>
 
 namespace qpragma::hhl::utils {

    /* Get the sign of a double*/
    double sign(const double & /*d*/);

    /* Normalize an array of coefficients */
    template <uint64_t SIZE>
    std::array<double, (1 << SIZE)> normalize(std::array<double, (1 << SIZE)> /*coeffs*/);

    /* Convert the binary value to a double flotting point value */
    double bin_to_double(uint64_t /*nb_bits*/, uint64_t /*val*/);

    /* Get the mean of each bit of all the eigenvalues */
    template <uint64_t SIZEC>
    std::array<double, SIZEC> get_means(std::vector<uint64_t> /*eigenvals*/);
 
 }   // qpragma::hhl::utils
 
 #include "utils.ipp"
 
 #endif  /* QPRAGMA_HHL_UTILS_HPP */
 