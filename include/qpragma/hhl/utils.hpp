/* -*- coding: utf-8 -*- */
/*
 * @file        qpragma/hhl/utils.hpp
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
