#include <cmath>
#include <unordered_map>
#include <complex>
#include <iostream>
// Include Q-Pragma framework
#include "qpragma.h"

// Include cpp files
// #include "utils.cpp"
// #include "observable.cpp"
#include "hamiltonian_sim.cpp"
#include "state_prep.cpp"


// Use namespace "qpragma"
using namespace qpragma;


const uint64_t NB_SHOTS = 50UL;
const double EPS = 1e-2;
const uint64_t NMAX = 100UL;
const std::string HAM_SIM_METHOD = "trotter";


/* Controlled UA part of the QPE */
#pragma quantum routine (Observable<SIZE> obs)
template <uint64_t SIZE, uint64_t SIZEC>
void ControlledUA(const quint_t<SIZE> & breg, const quint_t<SIZEC> & creg) {
    for (uint64_t i = 0 ; i < SIZEC ; ++i) {
        #pragma quantum ctrl (creg[SIZEC - i - 1])
        hamiltonian_sim<SIZE>(obs, -1, EPS, i, HAM_SIM_METHOD)(breg);
    }
}

/* Quantum Phase Estimation algorithm */
#pragma quantum routine (Observable<SIZE> obs)
template <uint64_t SIZE, uint64_t SIZEC>
void QPE(const quint_t<SIZE> & breg, const quint_t<SIZEC> & creg) {
    wall::H<SIZEC>(creg);
    (ControlledUA<SIZE, SIZEC>(obs))(breg, creg);
    qft<SIZEC>(creg);
}

/* Quantum Phase Estimation with 1 measurement of eigenvalue */
template <uint64_t SIZE, uint64_t SIZEC>
uint64_t QPEA(Observable<SIZE> obs, const quint_t<SIZE> & breg) {
    quint_t<SIZEC> creg;
    (QPE<SIZE, SIZEC>(obs))(breg, creg);
    return measure_and_reset(creg);
}

/* Get an estimation of all the eigenvalues by sampling on QPE */
template <uint64_t SIZE, uint64_t SIZEC>
std::vector<uint64_t> get_eigenvals(Observable<SIZE> obs,
                                    std::array<double, (1 << SIZE)> coeffs_b,
                                    uint64_t nb_shots=NB_SHOTS) {
    quint_t<SIZE> breg;
    // Initialize an array to check if a value is eigenvalue
    std::array<bool, 1 << SIZEC> is_eigen;
    for (int i = 0 ; i < (1 << SIZEC) ; ++i) {
        is_eigen[i] = false;
    }
    // Sample on QPE
    for (int i = 0 ; i < nb_shots ; ++i) {
        (state_prep<SIZE>(breg, coeffs_b));
        uint64_t res = QPEA<SIZE, SIZEC>(obs, breg);
        is_eigen[res] = true;
        reset(breg);
    }
    // Create the vector containing the eigenvalues
    std::vector<uint64_t> eigenvals;
    for (int i = 0 ; i < (1 << SIZEC) ; ++i) {
        if (is_eigen[i]) {
            eigenvals.push_back(i);
        }
    }

    return eigenvals;
}

/* Check if a value lambda is compatible with the means observed on eigenvalues */
template <uint64_t SIZEC>
bool is_compatible(uint64_t lambda, std::array<double, SIZEC> means, const uint64_t size) {
    for (int i = 0 ; i < size ; ++i) {
        if (means[i] == 0. or means[i] == 1.) {
            if (((lambda >> i) & 1) != means[i]) {
                return false;  // Is not compatible with the eigenvalues
            }
        }
    }
    return true;  // Is compatible
}

/* Reduced version of the AQE */
// Maybe this part does not work properly
#pragma quantum routine (double c, std::array<double, SIZEC> means)
template <uint64_t SIZEC>
void reduced_AQE(const quint_t<SIZEC> & creg, const qbool & anc) {
    for (uint64_t val_c = 0 ; val_c < (1 << SIZEC) ; ++val_c) {
        // Only on the compatible with eigenvalues
        if (is_compatible<SIZEC>(val_c, means, SIZEC)) {
            // Angle of the rotation RY
            double theta = 2. * acos(sqrt(1 - c*c /(val_c*val_c)));
            // Rotation on the ancilla controlled by the eigenvalue
            #pragma quantum ctrl (creg == val_c)
            (RY(theta))(anc);
        }
    }
}

/* Reduced version of the QFT */
#pragma quantum routine (std::array<double, SIZEC> means)
template <uint64_t SIZEC>
void reduced_qft(const quint_t<SIZEC> & creg) {
    for (int target = 0 ; target < SIZEC ; ++target) {
        // On 0 nothing happens here
        if (means[target] != 0.) {
            // Apply an H gate where the mean is not fixed
            if (means[target] != 1.) {
                H(creg[target]);
            }

            // controlled phase part of the QFT
            for (int control = target + 1 ; control < SIZEC ; ++control) {
                double angle = M_PI / (1 << (control - target));
                if (means[control] == 1.) {
                    // No control needed because control qubit always 1
                    (PH(angle))(creg[target]);
                }

                else if (means[control] != 0.) {
                    #pragma quantum ctrl (creg[control])
                    (PH(angle))(creg[target]);
                }
            }
        }
    }
}

/* Reduced QPE */
#pragma quantum routine (std::array<double, SIZEC> means, Observable<SIZE> obs)
template <uint64_t SIZE, uint64_t SIZEC>
void reduced_QPE(const quint_t<SIZE> & breg, const quint_t<SIZEC> & creg) {
   for (int i = 0 ; i < SIZEC ; ++i) {
       uint64_t idx = SIZEC - i - 1;
        // if qubit i in 1 --> X
        if (means[idx] == 1.) {
            X(creg[idx]);
            #pragma quantum ctrl (creg[idx])
            hamiltonian_sim<SIZE>(obs, -1, EPS, i, HAM_SIM_METHOD)(breg);
        }
        // if qubit mean not fixed --> H
        else if (means[i] != 0.) {
            H(creg[idx]);
            #pragma quantum ctrl (creg[idx])
            hamiltonian_sim<SIZE>(obs, -1, EPS, i, HAM_SIM_METHOD)(breg);
        }
        // if qubit i in 0 --> nothing happens
    }
    // Call the reduced QFT
    (reduced_qft<SIZEC>(means))(creg);
}

/* Reduced version of HHL */
template <uint64_t SIZE, uint64_t SIZEC>
void reduced_HHL(Observable<SIZE> obs, std::array<double, (1 << SIZE)> coeffs_b,
                 std::vector<uint64_t> eigenvals, double c,
                 quint_t<SIZE> & breg, quint_t<SIZEC> & creg) {
    qbool anc;
    // Get the means of the bits of the eigenvalues
    std::array<double, SIZEC> means = get_means<SIZEC>(eigenvals);
    // Post selection on ancilla in 1 state
    do {
        reset(breg);
        (state_prep<SIZE>(breg, coeffs_b))(breg);
        {
            #pragma quantum compute
            {
            (reduced_QPE<SIZE, SIZEC>(means, obs))(breg, creg);
            }
            
            (reduced_AQE<SIZEC>(c, means))(creg, anc);
        }
        // Automatically uncompute reduced_QPE
        
        reset(creg);
    } while (not measure_and_reset(anc));
}

/* Hybrid HHL algorithm */
template <uint64_t SIZE, uint64_t SIZEC>
void hybrid_HHL(Observable<SIZE> obs, std::array<double, (1 << SIZE)> coeffs_b,
                double c, const quint_t<SIZEC> & creg, const quint_t<SIZE> & breg,
                std::string ham_sim_method) {
    // Get the estimation of the eigenvalues of the matrix
    std::vector<uint64_t> eigenvals = get_eigenvals<SIZE, SIZEC>(obs, coeffs_b);
    // Call reduced HHL
    reduced_HHL<SIZE, SIZEC>(obs, coeffs_b, eigenvals, c, breg, creg);
}
