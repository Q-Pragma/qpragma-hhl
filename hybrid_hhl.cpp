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
const uint64_t N_ITER_TROTTER = 100UL;
const uint64_t NMAX = 100UL;

double bin_to_double(uint64_t nb_bits, uint64_t val) {
    if (val == 0) {
        return 1.;
    }
    else {
        uint64_t reverse = 0;
        for (int i = nb_bits - 1 ; i >=0 ; --i) {
            reverse |= (val & 1UL) << i;
            val >>= 1;
        }
        return (double) reverse / (double) (1 << nb_bits);
    }
}

#pragma quantum routine (std::array<double, 2UL> coeffs)
void state_prep(const qbool & q) {
    // Prepare the superposition of coeffs a and b in qreg
    // a|0> + b|1>
    double a = coeffs[0];
    RY(2 * acos(a))(q);
}

/* Hamiltonian simulation of A using Trotterization
 * UA = e^(2i*pi*A) + O(...)
*/
#pragma quantum routine (uint64_t n, std::array<double, 4UL> coeffs)
void UA (const qbool & breg) {
    // Identity --> just add a global phase
    qbool anc;
    CNOT(breg, anc);
    RZ(4. * M_PI * coeffs[3] * (1 << n)).ctrl(anc, breg);
    RZ(-4. * M_PI * coeffs[3] * (1 << n)).ctrl((qbool) not anc, breg);
    CNOT(breg,anc);
    // Special case: only X and I
    if (coeffs[1] == 0. and coeffs[2] == 0.) {
        RX(-4. * M_PI * coeffs[0] * (1 << n))(breg);
    }
    // Special case: only Y and I
    else if (coeffs[0] == 0. and coeffs[2] == 0.) {
        RY(-4. * M_PI * coeffs[1] * (1 << n))(breg);
    }
    // Special case: only Z and I
    else if (coeffs[0] == 0. and coeffs[1] == 0.) {
        RZ(-4. * M_PI * coeffs[2] * (1 << n))(breg);
    }
    else {
        for (int i = 0 ; i < (1 << n) ; ++i) {
            for (int j = 0 ; j < N_ITER_TROTTER ; ++j) {
                RX(-4. * M_PI * coeffs[0] / N_ITER_TROTTER)(breg);
                RY(-4. * M_PI * coeffs[1] / N_ITER_TROTTER)(breg);
                RZ(-4. * M_PI * coeffs[2] / N_ITER_TROTTER)(breg);
            }
        }
    }
}

# pragma quantum routine (std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void ControlledUA(const qbool & breg, const quint_t<SIZEC> & creg) {
    for (uint64_t i = 0 ; i < SIZEC ; ++i) {
        (UA(i, coeffs_mat)).ctrl(creg[i], breg);
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

            prev = curr;
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
    for (uint64_t val_c = 0 ; val_c < (1 << SIZEC) ; ++val_c) {
        if (is_compatible(val_c, fixed_idx)) {
            double lambda = bin_to_double(SIZEC, val_c);
            theta = c / lambda;
            theta = 2 * acos( sqrt(1 - theta * theta) );
            // Optimizable with greycode
            #pragma quantum ctrl (creg == val_c)
            (RY(theta))(anc);
        }
    }
}

#pragma quantum routine (map_idmean fixed_idx, std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void reduced_QPE(const qbool & breg, const quint_t<SIZEC> & creg) {
   wall::H<SIZEC>(creg);
   for (auto idx_m : fixed_idx) {
       uint64_t idx = idx_m.first;
       bool m = idx_m.second;
       if (m) {
           X(creg[idx]);
       }
       H(creg[idx]);
       (UA(idx, coeffs_mat)).ctrl(creg[idx], breg);
   }
   qft<SIZEC>.dag(creg);
}

template <uint64_t SIZEC>
void reduced_HHL(std::array<double, 4UL> coeffs_mat, distribution & distr, double c,
                 const qbool & breg, const quint_t<SIZEC> & creg) {
    qbool anc;
    map_idmean fixed_idx = get_fixed(distr, SIZEC);
    // 1 state is never reached :(
    //do {
        (reduced_QPE<SIZEC>(fixed_idx, coeffs_mat))(breg, creg);
        (reduced_AQE<SIZEC>(c, fixed_idx))(creg, anc);
        (reduced_QPE<SIZEC>(fixed_idx, coeffs_mat)).dag(breg, creg);
    //} while (not measure_and_reset(anc));
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

template <uint64_t SIZEC>
void bug_solver(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                 uint64_t nb_shots=NB_SHOTS) {
    qbool breg;
    quint_t<SIZEC> creg;
    uint64_t res[2UL] = {0, 0};
    for (int i = 0 ; i < nb_shots ; ++i) {
        hybrid_HHL<SIZEC>(coeffs_mat, coeffs_b, creg, breg);
        uint64_t valc = measure_and_reset(creg);
        H(breg);
        bool valb = measure_and_reset(breg);
        ++res[valb];
    }
    std::cout << "Final result: " << std::endl;
    std::cout << "0: " << res[0] / (double) nb_shots << " 1: "
              << res[1] / (double) nb_shots << std::endl;
}


template <uint64_t SIZEC>
void test_solver(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                 uint64_t nb_shots=NB_SHOTS) {
    double c = std::abs(condition_number(coeffs_mat));
    distribution distr = get_distribution<SIZEC>(coeffs_mat, coeffs_b);
    qbool breg;
    quint_t<SIZEC> creg;
    uint64_t res[2UL] = {0, 0};
    bool valb;
    uint64_t valc;
    for (int i = 0 ; i < nb_shots ; ++i) {
        do {
            (state_prep(coeffs_b))(breg);
            reduced_HHL<SIZEC>(coeffs_mat, distr, c ,breg, creg);
            valc = measure_and_reset(creg);
            H(breg);
            valb = measure_and_reset(breg);
        } while (valc != 0);
        ++res[valb];
    }
    std::cout << "Final result: " << std::endl;
    std::cout << "0: " << res[0] / (double) nb_shots << " 1: "
              << res[1] / (double) nb_shots << std::endl;
}

// Main function
int main() {
    const uint64_t SIZEC = 3UL;
    qbool breg;
    quint_t<SIZEC> creg;
    std::array<double, 4UL> coeffs_mat = {0.25,0.,0.,0.5};
    //coeffs_mat = normalize<4UL>(coeffs_mat);
    std::array<double, 2UL> coeffs_b = {1.,0.};
    coeffs_b = normalize<2UL>(coeffs_b);
    uint64_t NB_SHOT = 10UL;
    std::vector<uint64_t> res(1<<SIZEC);
    #pragma quantum scope
    // Display the eigenvalues of A
    // /!\ They have to be in }0,1)
    {
    for (int i = 0 ; i < 1000 ; ++i) {
        (QPE<SIZEC>(coeffs_mat))(breg, creg);
        res[measure_and_reset(creg)]++;
        reset(breg);
    }
    for (int i = 0 ; i < (1 << SIZEC) ; ++i) {
        std::cout << bin_to_double(SIZEC, i) << ":" << res[i] << std::endl;
    }
    // Test hhl: Do not work properly for the moment
    test_solver<SIZEC>(coeffs_mat, coeffs_b, NB_SHOT);
    }
}
