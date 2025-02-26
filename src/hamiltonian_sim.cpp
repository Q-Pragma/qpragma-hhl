#include <array>
#include <random>
#include <iostream>

#include "qpragma.h"
//#include "qsvt.cpp"
#include "observable.cpp"

using namespace qpragma;

const uint64_t N_ITER_TROTTER = 10UL;

/**
 * Hamiltonian simulation using QSVT
 */
/*
#pragma quantum routine (std::array<stad::array<double, N>, N> matrix, double t)
template <uint64_t N>
void hamiltonian_simulation_via_qsvt(const array<N> & target) {
    qbool anc = plus;
    qsvt::qsvt_cos(matrix, t).ctrl(not anc, target);
    qsvt::qsvt_sin(matrix, t).ctrl(anc, target);
}
*/




/**
 * Hamiltonian simulation via Trotterization
 */

 #pragma quantum routine (Observable<SIZE> obs, double t, double eps, uint64_t pow2)
 template <uint64_t SIZE>
 void hamiltonian_sim_via_trotter(const array<SIZE> & qreg) {
    uint64_t n_trotter = N_ITER_TROTTER;  // TODO
    // Create a new observable with all the coeffs * -t / n_trotter
    std::vector<std::string> strings = obs.get_strings();
    std::vector<double> coeffs = obs.get_coeffs();
    std::transform(coeffs.begin(), coeffs.end(),
                   coeffs.begin(),
                   [n_trotter, t](double c){ return  -t * c / (double) n_trotter; });
    Observable<SIZE> trotter_obs(strings, coeffs);

    // Trotterization
    for (uint64_t power = 0 ; power < (1 << pow2) ; ++power) {
        for (uint64_t step = 0 ; step < n_trotter ; ++step) {
            for (PauliTerm<SIZE> term : trotter_obs.pauli_terms) {
                (apply_term<SIZE>(term))(qreg);
            }
        }
    }
 }





/**
 * Hamiltonian simulation via qDrift
 */

 uint64_t qdrift_sampling(const std::vector<double> & weights) {
    // Get the distribution
    std::vector<double> distr(weights.size());
    double cum_weight;
    for (uint64_t idx = 0UL ; idx < weights.size() ; ++idx) {
        cum_weight += weights[idx];
        distr[idx] = cum_weight;
    }

    // Perform the sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0., cum_weight);
    double sample = dist(gen);

    for (uint64_t idx = 0UL ; idx < weights.size() ; ++idx) {
        if (distr[idx] >= sample) return idx;
    }
    return weights.size();
}

template <uint64_t SIZE>
std::vector<PauliTerm<SIZE>> qdrift_unitary_list(Observable<SIZE> obs,
                                                 double t, double eps) {
    std::vector<std::string> strings = obs.get_strings();
    std::vector<double> coeffs = obs.get_coeffs();

    double lambda = std::accumulate(coeffs.cbegin(), coeffs.cend(), 0.);
    uint64_t N = ceil(2. * lambda * lambda * t * t / eps);

    std::vector<PauliTerm<SIZE>> unitaries;
    for (uint64_t i = 0UL ; i < N ; ++i) {
        uint64_t idx = qdrift_sampling(coeffs);
        unitaries.push_back(PauliTerm<SIZE>(strings[idx], lambda * t / N));
    }

    return unitaries;
}

#pragma quantum routine (Observable<SIZE> obs, double t, double eps, uint64_t pow2)
template <uint64_t SIZE>
void hamiltonian_sim_via_qdrift(const array<SIZE> & qreg) {
    std::vector<PauliTerm<SIZE>> unitaries = qdrift_unitary_list(obs, t, eps);
    for (uint64_t power = 0 ; power < (1 << pow2) ; ++power) {
        for (int i = 0 ; i < obs.pauli_terms.size() ; ++i) {
            apply_term<SIZE>(obs.pauli_terms[i])(qreg);
        }
    }
}





/* Hamiltonian simulation of A^{2^n} using Trotterization
 * UA = e^(2i*pi*A) + O(...)
*/
/*
#pragma quantum routine (const observable::Observable & H, uint64_t pow2)
void hamiltonian_sim_2x2(const qbool & breg) {
    std::vector<double> coeffs = H.get_coeffs();
    std::vector<std::string> strings = H.get_string();


    // Identity --> just add a global phase
    qbool anc;
    CNOT(breg, anc);
    RZ(-4. * M_PI * coeffs[3] * (1 << n)).ctrl(anc, breg);
    RZ(4. * M_PI * coeffs[3] * (1 << n)).ctrl((qbool) not anc, breg);
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
    // General case
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
    */


#pragma quantum routine (Observable<SIZE> obs, double t, double eps, uint64_t pow2, std::string method)
template <uint64_t SIZE>
void hamiltonian_sim(const array<SIZE> & qreg) {
    if (method == "qdrift") {
        hamiltonian_sim_via_qdrift<SIZE>(obs, t, eps, pow2)(qreg);
    }
    else if (method == "trotter") {
        hamiltonian_sim_via_trotter<SIZE>(obs, t, eps, pow2)(qreg);
    }
}