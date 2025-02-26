#include <vector>
#include <string>
#include <complex>
#include <cstdint>

#include "qpragma.h"

using namespace qpragma;


template <uint64_t SIZE>
class PauliTerm {
public:
    std::string pauli_string;
    double coeff;

    PauliTerm(std::string pauli_string, double coeff)
        : pauli_string(pauli_string), coeff(coeff) {
            // TODO ASSERT on SIZE
        }

    PauliTerm(const PauliTerm<SIZE> & other) 
        : pauli_string(other.pauli_string), coeff(other.coeff) {}
};


template <uint64_t SIZE>
class Observable {
public:
    std::vector<PauliTerm<SIZE>> pauli_terms;

    Observable(std::vector<PauliTerm<SIZE>> pauli_terms)
        : pauli_terms(pauli_terms) {}

    Observable(std::vector<std::string> pauli_strings, std::vector<double> coeffs) {
            // TODO ASSERT on SIZE
            for (uint64_t idx = 0UL ; idx < coeffs.size() ; ++idx) {
                pauli_terms.push_back(PauliTerm<SIZE>(pauli_strings[idx], coeffs[idx]));
            }
    }

    std::vector<double> get_coeffs() const {
        std::vector<double> coeffs(pauli_terms.size());
        std::transform(pauli_terms.begin(), pauli_terms.end(),
                        coeffs.begin(),
                        [](PauliTerm<SIZE> p){ return p.coeff; });
        return coeffs;
    }

    std::vector<std::string> get_strings() const {
        std::vector<std::string> strings;
        for (uint64_t i = 0 ; i < pauli_terms.size() ; ++i) 
            strings.push_back(pauli_terms[i].pauli_string);
        return strings;
    }
};

#pragma quantum routine (PauliTerm<SIZE> term)
template <uint64_t SIZE>
void apply_term(const array<SIZE> & qreg) {
    if (std::count(term.pauli_string.begin() , term.pauli_string.end(), 'I') == SIZE) {
        // TODO : Add a global phase
    } else {
        qbool anc;
        # pragma quantum compute
        {
            for (uint64_t idx = 0UL ; idx < SIZE ; ++idx) {
                // I --> Do nothing
                if (term.pauli_string[idx] == 'I') {
                    continue;
                }
                // X = H Z H
                if (term.pauli_string[idx] == 'X') {
                    H(qreg[idx]);
                }
                // Y = Sdag H Z H S 
                else if (term.pauli_string[idx] == 'Y') {
                    RX(-M_PI / 2.)(qreg[idx]);
                }
                // Stock the parity on XYZ qubits in anc
                CNOT(qreg[idx], anc);
            }
        }
        RZ(term.coeff)(anc);
    }
}
