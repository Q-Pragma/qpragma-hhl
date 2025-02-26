#include "hybrid_hhl.cpp"


/* Test the solver implementation */
template <uint64_t SIZE, uint64_t SIZEC>
void test_solver(Observable<SIZE> obs, std::array<double, (1 << SIZE)> coeffs_b,
                 double c, uint64_t nb_shots=NB_SHOTS) {
    // Get the eigenvalues of the matrix
    std::vector<uint64_t> eigenvals = get_eigenvals<SIZE, SIZEC>(obs, coeffs_b);
    // Print them
    std::cout << "Eigenvalues" << std::endl;
    for (auto l : eigenvals) {
        std::cout << "binary: " << l
                  << ", decimal: " << bin_to_double(SIZEC, l)
                  << std::endl;
    }
    std::cout << "-----------------------" << std::endl;
    /*
    // Allocate quantum registers
    quint_t<SIZE> breg;
    quint_t<SIZEC> creg;

    // Final result probabilities
    std::vector<uint64_t> res(1 << SIZE);


    uint64_t valb;
    uint64_t valc;
    for (int i = 0 ; i < nb_shots ; ++i) {
        // Call reduced HHL
        reduced_HHL<SIZE, SIZEC>(H, coeffs_b, eigenvals, c, breg, creg);
        // Measure creg in Z basis
        valc = measure_and_reset(creg);
        // Measure breg in X basis
        valb = measure_and_reset(breg);
        // Store the measurement result
        ++res[valb];
    }

    std::cout << "Final result: " << std::endl;
    // TODO
    */
}

// Main function
int main() {
    const uint64_t SIZEC = 2UL;
    const uint64_t SIZE = 1UL;
    quint_t<SIZE> breg;
    quint_t<SIZEC> creg;
    std::vector<std::string> pauli_strings({"X"});
    std::vector<double> coeffs({0.5});
    Observable<SIZE> obs(pauli_strings, coeffs);
    double c = 0.25/0.75;
    std::array<double, (1 << SIZE)> coeffs_b;
    coeffs_b[0] = 1.;
    coeffs_b = normalize<SIZE>(coeffs_b);
    uint64_t NB_SHOT = 1000UL;
    // Test hhl
    test_solver<SIZE, SIZEC>(obs, coeffs_b, c, NB_SHOT);
    //hamiltonian_sim_via_trotter<SIZE>(obs, -1, EPS, 0)(breg);
}
