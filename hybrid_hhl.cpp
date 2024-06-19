#include <cmath>
#include <unordered_map>
#include <complex>
#include <iostream>
// Include Q-Pragma framework
#include "qpragma.h"

// Use namespace "qpragma"
using namespace qpragma;


const uint64_t NB_SHOTS = 100UL;
const double EPS = 1e-3;
const uint64_t N_ITER_TROTTER = 10UL;
const uint64_t NMAX = 100UL;


/* Convert the binary value to a double flotting point value */
double bin_to_double(uint64_t nb_bits, uint64_t val) {
    if (val == 0) {
        return 1.;  // 0.000 and 1.000 equivalent but only 1. can be eigenvalue
    }
    else {
        return (double) val / (double) (1 << nb_bits);
    }
}

// Prepare the superposition of coeffs a and b in qreg
// a|0> + b|1>
#pragma quantum routine (std::array<double, 2UL> coeffs)
void state_prep(const qbool & q) {
    double a = coeffs[0];
    RY(2 * acos(a))(q);
}

/* Hamiltonian simulation of A^{2^n} using Trotterization
 * UA = e^(2i*pi*A) + O(...)
*/
#pragma quantum routine (uint64_t n, std::array<double, 4UL> coeffs)
void UA(const qbool & breg) {
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

/* Controlled UA part of the QPE */
#pragma quantum routine (std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void ControlledUA(const qbool & breg, const quint_t<SIZEC> & creg) {
    for (uint64_t i = 0 ; i < SIZEC ; ++i) {
        #pragma quantum ctrl (creg[SIZEC - i - 1])
        (UA(i, coeffs_mat))(breg);
    }
}

/* Quantum Phase Estimation algorithm */
#pragma quantum routine (std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void QPE(const qbool & breg, const quint_t<SIZEC> & creg) {
    wall::H<SIZEC>(creg);
    (ControlledUA<SIZEC>(coeffs_mat))(breg, creg);
    qft<SIZEC>(creg);
}

/* Quantum Phase Estimation with 1 measurement of eigenvalue */
template <uint64_t SIZEC>
uint64_t QPEA(std::array<double, 4UL> coeffs_mat, const qbool & breg) {
    quint_t<SIZEC> creg;
    (QPE<SIZEC>(coeffs_mat))(breg, creg);
    return measure_and_reset(creg);
}

/* Get an estimation of all the eigenvalues by sampling on QPE */
template <uint64_t SIZEC>
std::vector<uint64_t> get_eigenvals(std::array<double, 4UL> coeffs_mat,
                                    std::array<double, 2UL> coeffs_b,
                                    uint64_t nb_shots=NB_SHOTS) {
    qbool breg;
    // Initialize an array to check if a value is eigenvalue
    std::array<bool, 1 << SIZEC> is_eigen;
    for (int i = 0 ; i < (1 << SIZEC) ; ++i) {
        is_eigen[i] = false;
    }

    // Sample on QPE
    for (int i = 0 ; i < nb_shots ; ++i) {
        (state_prep(coeffs_b))(breg);
        uint64_t res = QPEA<SIZEC>(coeffs_mat, breg);
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

/* Get the mean of each bit of all the eigenvalues */
template <uint64_t SIZEC>
std::array<double, SIZEC> get_means(std::vector<uint64_t> eigenvals) {
    std::array<double, SIZEC> means;
    for (int i = 0 ; i < SIZEC ; ++i) {
        means[i] = 0.;
        // Compute the mean on bit i
        for (uint64_t val : eigenvals) {
            means[i] += (double) ((val >> i) & 1);
        }
        means[i] /= (double) eigenvals.size();
    }
    return means;
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
#pragma quantum routine (std::array<double, SIZEC> means, std::array<double, 4UL> coeffs_mat)
template <uint64_t SIZEC>
void reduced_QPE(const qbool & breg, const quint_t<SIZEC> & creg) {
   for (int i = 0 ; i < SIZEC ; ++i) {
       uint64_t idx = SIZEC - i - 1;
        // if qubit i in 1 --> X
        if (means[idx] == 1.) {
            X(creg[idx]);
            #pragma quantum ctrl (creg[idx])
            (UA(i, coeffs_mat))(breg);
        }
        // if qubit mean not fixed --> H
        else if (means[i] != 0.) {
            H(creg[idx]);
            #pragma quantum ctrl (creg[idx])
            (UA(i, coeffs_mat))(breg);
        }
        // if qubit i in 0 --> nothing happens
    }
    // Call the reduced QFT
    (reduced_qft<SIZEC>(means))(creg);
}

/* Reduced version of HHL */
template <uint64_t SIZEC>
void reduced_HHL(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                 std::vector<uint64_t> eigenvals, double c,
                 qbool & breg, quint_t<SIZEC> & creg) {
    qbool anc;
    // Get the means of the bits of the eigenvalues
    std::array<double, SIZEC> means = get_means<SIZEC>(eigenvals);
    // Post selection on ancilla in 1 state
    do {
        reset(breg);
        (state_prep(coeffs_b))(breg);
        {
        #pragma quantum compute
        (reduced_QPE<SIZEC>(means, coeffs_mat))(breg, creg);
        
        (reduced_AQE<SIZEC>(c, means))(creg, anc);
        }
        // Automatically uncompute reduced_QPE
        
        reset(creg);
    } while (not measure_and_reset(anc));
}

/* Normalize an array of coefficients */
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

/* Power method to find the highest eigenvalue */
std::complex<double> power_method(std::array<double, 4UL> coeffs_mat,
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

/* Hybrid HHL algorithm */
template <uint64_t SIZEC>
void hybrid_HHL(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                double c, const quint_t<SIZEC> & creg, const qbool & breg) {
    // Get the estimation of the eigenvalues of the matrix
    std::vector<uint64_t> eigenvals = get_eigenvals<SIZEC>(coeffs_mat, coeffs_b);
    // Call reduced HHL
    reduced_HHL<SIZEC>(coeffs_mat, coeffs_b, eigenvals, c, breg, creg);
}

/* Test the solver implementation */
template <uint64_t SIZEC>
void test_solver(std::array<double, 4UL> coeffs_mat, std::array<double, 2UL> coeffs_b,
                 double c, uint64_t nb_shots=NB_SHOTS) {
    // Get the eigenvalues of the matrix
    std::vector<uint64_t> eigenvals = get_eigenvals<SIZEC>(coeffs_mat, coeffs_b);
    // Print them
    std::cout << "Eigenvalues" << std::endl;
    for (auto l : eigenvals) {
        std::cout << "binary: " << l
                  << ", decimal: " << bin_to_double(SIZEC, l)
                  << std::endl;
    }
    std::cout << std::endl;

    // Allocate quantum registers
    qbool breg;
    quint_t<SIZEC> creg;

    // Final result probabilities
    uint64_t res[2UL] = {0, 0};


    bool valb;
    uint64_t valc;
    for (int i = 0 ; i < nb_shots ; ++i) {
        // Call reduced HHL
        reduced_HHL<SIZEC>(coeffs_mat, coeffs_b, eigenvals, c ,breg, creg);
        // Measure creg in Z basis
        valc = measure_and_reset(creg);
        // Measure breg in X basis
        valb = measure_and_reset(breg);
        // Store the measurement result
        ++res[valb];
    }

    std::cout << "Final result: " << std::endl;
    std::cout << "0: " << res[0] / (double) nb_shots << " 1: "
              << res[1] / (double) nb_shots << std::endl;
}

// Main function
int main() {
    const uint64_t SIZEC = 2UL;
    qbool breg;
    quint_t<SIZEC> creg;
    std::array<double, 4UL> coeffs_mat = {0.25,0.,0.,0.5};
    double c = 0.25/0.75;
    std::array<double, 2UL> coeffs_b = {1., 0.};
    coeffs_b = normalize<2UL>(coeffs_b);
    uint64_t NB_SHOT = 1000UL;
    std::vector<uint64_t> res(1<<SIZEC);
    #pragma quantum scope
    {
    // Test hhl
    test_solver<SIZEC>(coeffs_mat, coeffs_b, c, NB_SHOT);
    }
}
