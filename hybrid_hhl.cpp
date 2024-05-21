#include <cmath>
#include <unordered_map>
#include <complex>
#include <iostream>
// Include Q-Pragma framework
#include "qpragma.h"

// Use namespace "qpragma"
using namespace qpragma;


using distribution = std::unordered_map<uint64_t, double>;
using map_idmean = std::unordered_map<uint64_t, bool>;

const uint64_t NB_SHOTS = 1000UL;
const double EPS = 1e-3;
const double NMAX = 100000;


#pragma quantum routine (std::array<double, 2UL> coeffs)
void state_prep(const qbool & q) {
    // Prepare the superposition of coeffs a and b in qreg
    // a|0> + b|1>
    double a = coeffs[0];
    RY(2 * acos(a))(q);
}

#pragma quantum routine (uint64_t n, std::array<double, 4UL> coeffs)
void UA(const qbool & breg) {
    // Special cases : already a Pauli
    // Pauli X
    if (coeffs[0] == 1.) {if (n == 0) {X(breg);}
        return;
    }
    // Pauli Y
    else if (coeffs[1] == 1.) {
        if (n == 0) {Y(breg);}
        return;
    }
    // Pauli Z
    else if (coeffs[2] == 1.) {
        if (n == 0) {Z(breg);}
        return;
    }
    else if (coeffs[3] == 1.) {
        return;
    }


    // General case
    // The operator UA can be decomposed using rotations
    // UA = e^(i alpha) RZ(theta2) RY(theta1) RZ(theta0)
    std::complex<double> U[2][2] = {
        {coeffs[2] + coeffs[3], coeffs[0] + 1i * coeffs[1]},
        {coeffs[0] + 1i * coeffs[1], coeffs[3] - coeffs[2]}
    };
    // Compute alpha
    std::complex<double> detU = U[0][0] * U[1][1] - U[0][1] * U[1][0];
    double alpha = 0.5 * atan2(std::imag(detU), std::real(detU));
    for (int i = 0 ; i < 2UL ; ++i){
        for (int j = 0 ; j < 2UL ; ++j) {
            std::complex<double> factor = cos(alpha) + 1i * sin(alpha);
            U[i][j] *= factor;
        }                                                                         
    }
    // Compute theta0, theta1 and theta2
    double theta0, theta1, theta2;
    if (std::abs(U[0][0]) > std::abs(U[0][1])) {
        theta1 = 2. * acos(std::abs(U[0][0]));
    }
    else {
        theta1 = 2. * asin(std::abs(U[0][1]));
    }
    double sum = 0.;
    if (cos(theta1 / 2.) != 0.) {
        std::complex<double> val = U[1][1] / cos(theta1 / 2.);
        sum = 2 * atan2(std::imag(val), std::real(val));
    }
    double diff = 0.;
    if (sin(theta1 / 2.) != 0.) {
        std::complex<double> val = U[1][0] / sin(theta1 / 2.);
        diff = 2 * atan2(std::imag(val), std::real(val));
    }
    theta0 = (sum + diff) / 2.;
    theta2 = (sum - diff) / 2.;


    // Build the corresponding quantum circuit
    (RZ(theta0))(breg);  // First RZ(theta0) not merged with RZ(theta2)
    (RY(theta1))(breg);  // First RY out of for loop
    for (int i = 1 ; i < (1 << n) ; ++i) {
        (RZ(theta2 + theta0))(breg);  // Merge theta0 and theta2
        (RY(theta1))(breg);
    }
    (RZ(theta2))(breg);  // Last RZ(theta2) not merged with RZ(theta0)
    PH(std::pow(alpha, (1 << n)))(breg);  // PH commutes with everything before
}

# pragma quantum routine (std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void ControlledUA(const qbool & breg, const quint_t<SIZEC> & creg) {
    for (uint64_t i = 0 ; i < SIZEC ; ++i) {
        (UA(i, coeffs_mat)).ctrl(creg[SIZEC-i-1], breg);
    }
}

#pragma quantum routine (std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void QPE(const qbool & breg, const quint_t<SIZEC> & creg) {
    wall::H<SIZEC>(creg);
    (ControlledUA<SIZEC>(coeffs_mat))(breg, creg);
    qft<SIZEC>.dag(creg);
}

template <uint64_t SIZEC>
uint64_t QPEA(std::array<double, 4UL> coeffs_mat, const qbool & breg) {
    quint_t<SIZEC> creg;
    (QPE<SIZEC>(coeffs_mat))(breg, creg);
    return measure_and_reset(creg);
}

template <uint64_t SIZEC>
distribution get_distribution(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                              uint64_t nb_shots=NB_SHOTS) {
    qbool breg;
    distribution distr;
    double increment = 1. / (double) nb_shots;
    for (int i = 0 ; i < nb_shots ; ++i) {
        (state_prep(coeffs_b))(breg);
        uint64_t res = QPEA<SIZEC>(coeffs_mat, breg);
        distr[res] += increment;
        reset(breg);
    }
    return distr;
}

map_idmean get_fixed(distribution distr, uint64_t size) {
    map_idmean fixed_idx;
    for (uint64_t idx = 0UL ; idx < size ; ++idx) {
        int prev = -1;
        int curr;
        bool fixed = true;

        for (auto d : distr) {
            curr = d.first & (1 << idx);

            if (prev != -1) {  
                fixed = !(prev ^ curr);
            }

            if (!fixed){
                break;
            }
        }

        if (fixed) {
            fixed_idx[idx] = curr;
        }
    }
    return fixed_idx;
}

bool is_compatible(uint64_t lambda, map_idmean fixed_idx) {
    for (auto e : fixed_idx) {
        if (not (((lambda >> e.first) & 1UL) == e.second) ) {
            return false;
        }
    }
    return true;
}

#pragma quantum routine (double c, map_idmean fixed_idx)
template <uint64_t SIZEC>
void reduced_AQE(const quint_t<SIZEC> & creg, const qbool & anc) {
    double theta;
    for (uint64_t lambda = 0 ; lambda < (1 << SIZEC) ; ++lambda) {
        if (is_compatible(lambda, fixed_idx)) {
            theta = c / (double) lambda;
            theta = 2 * acos( sqrt(1 - theta * theta) );
            // Optimizable with greycode
            #pragma quantum ctrl (creg == lambda)
            (RY(theta))(anc);
        }
    }
}


#pragma quantum routine (map_idmean fixed_idx, std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void reduced_QPE(const qbool & breg, const quint_t<SIZEC> & creg) {
    for (auto idx_m : fixed_idx) {  
        uint64_t idx = idx_m.first;
        bool m = idx_m.second;
        if (m) {
            X(creg[idx]);
        }
        H(creg[idx]);
        (UA(idx, coeffs_mat)).ctrl(creg[SIZEC - idx - 1], breg);
    }
   qft<SIZEC>.dag(creg);
}

template <uint64_t SIZEC>
void reduced_HHL(std::array<double, 4UL> coeffs_mat, distribution & distr, double c,
                 const qbool & breg, const quint_t<SIZEC> & creg) {
    qbool anc;
    map_idmean fixed_idx = get_fixed(distr, SIZEC);
    (reduced_QPE<SIZEC>(fixed_idx, coeffs_mat))(breg, creg);
    (reduced_AQE<SIZEC>(c, fixed_idx))(creg, anc);
    (reduced_QPE<SIZEC>(fixed_idx, coeffs_mat)).dag(breg, creg);
    reset(anc);
}

template <uint64_t SIZE>
std::array<double, SIZE> normalize(std::array<double, SIZE> coeffs) {
    double norm = 0.;
    std::array<double, SIZE> res;
    for (int i = 0 ; i < SIZE ; ++i) {
        norm += coeffs[i] * coeffs[i];
    }
    norm = sqrt(norm);
    for (int i = 0 ; i < SIZE ; ++i) {
        res[i] = coeffs[i] / norm;
    }
    return res;
}

std::complex<double> condition_number(std::array<double, 4UL> coeffs_mat,
                                      double n_iter=NMAX) {
    std::complex<double> q[2] = {1., 1.};
    std::complex<double> lambda;
    std::complex<double> z[2];
    for (int i = 0 ; i < n_iter ; ++i) {
        // Compute z
        z[0] = q[0] * std::complex<double>(coeffs_mat[3] + coeffs_mat[2])
             + q[1] * std::complex<double>(coeffs_mat[0] - 1i * coeffs_mat[1]);
        z[1] = q[0] * std::complex<double>(coeffs_mat[0] + 1i * coeffs_mat[1])
             + q[1] * std::complex<double>(coeffs_mat[3] - coeffs_mat[2]);  
        // Compute lambda
        lambda = q[0] * z[0] + q[1] * z[1];
        // Compute q
        double norm = sqrt( std::norm(z[0]) + std::norm(z[1]) );
        q[0] = z[0] / norm;
        q[1] = z[1] / norm;
    }
    return lambda;
}


template <uint64_t SIZEC>
void hybrid_HHL(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                const quint_t<SIZEC> & creg, const qbool & breg) {
    double c = std::abs(condition_number(coeffs_mat));
    distribution distr = get_distribution<SIZEC>(coeffs_mat, coeffs_b);
    // Display the distribution obtained with QPEA
    // Problem: after a first call to hybrid_HHL, the distribution got is always
    // broken (just 0 is obtained) --> I cannot explain it for the moment
    for (auto e : distr) {
        std::cout << e.first << " " << e.second << std::endl;
    }
    std::cout << "------" << std::endl;
    (state_prep(coeffs_b))(breg);
    reduced_HHL<SIZEC>(coeffs_mat, distr, c, breg, creg);
}


// Main function
int main() {
    const uint64_t SIZEC = 2UL;
    qbool breg;
    quint_t<SIZEC> creg;
    std::array<double, 4UL> coeffs_mat = {0.25,0.,0.,0.5};
    coeffs_mat = normalize<4UL>(coeffs_mat);
    std::array<double, 2UL> coeffs_b = {1.,0.};
    uint64_t NB_SHOT = 10UL;
    uint64_t res[2UL] = {0,0};
    for (int i = 0; i < NB_SHOT ; ++i) {
        hybrid_HHL<SIZEC>(coeffs_mat, coeffs_b, creg, breg);
        auto valc = measure_and_reset(creg);
        X(breg);
        bool valb = measure_and_reset(breg);
        ++res[valb];
    }
    std::cout << "Final result:" << std::endl;
    std::cout << "0: " << res[0] / (float) NB_SHOT << " 1: " << res[1] / (float) NB_SHOT << std::endl;
}
