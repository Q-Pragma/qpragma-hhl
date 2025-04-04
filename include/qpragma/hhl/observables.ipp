// #include "qpragma/hhl/observables.hpp"


// Pauli term
template <uint64_t SIZE>
inline qpragma::hhl::observables::PauliTerm<SIZE>::PauliTerm(
    std::span<const char> pauli_string, double coeff
): m_coeff(coeff) {
    std::transform(pauli_string.begin(), pauli_string.begin() + SIZE, m_pauli_string.begin(),
                   [](char op){ return static_cast<pauli_op>(op); });
}


template <uint64_t SIZE>
inline qpragma::hhl::observables::PauliTerm<SIZE>::PauliTerm(
    const std::array<qpragma::hhl::observables::pauli_op, SIZE> & pauli_string, double coeff
): m_pauli_string(pauli_string), m_coeff(coeff) {}



template <uint64_t SIZE>
inline double qpragma::hhl::observables::PauliTerm<SIZE>::coeff() const {
    return m_coeff;
}


template <uint64_t SIZE>
inline std::array<qpragma::hhl::observables::pauli_op, SIZE>
qpragma::hhl::observables::PauliTerm<SIZE>::pauli_string() const {
    return m_pauli_string;
}


template <uint64_t SIZE>
inline std::array<qpragma::hhl::observables::pauli_op, SIZE>::const_iterator
qpragma::hhl::observables::PauliTerm<SIZE>::begin() const {
    return m_pauli_string.cbegin();
}


template <uint64_t SIZE>
inline std::array<qpragma::hhl::observables::pauli_op, SIZE>::const_iterator
qpragma::hhl::observables::PauliTerm<SIZE>::end() const {
    return m_pauli_string.cend();
}


template <uint64_t SIZE>
inline bool qpragma::hhl::observables::PauliTerm<SIZE>::is_identity() const {
    return std::count(m_pauli_string.begin(), m_pauli_string.end(), pauli_op::I) == SIZE;
}


template <uint64_t SIZE>
inline qpragma::hhl::observables::pauli_op qpragma::hhl::observables::PauliTerm<SIZE>::operator[](uint64_t idx) const {
    return m_pauli_string[idx];
}


// Observable
template <uint64_t SIZE>
inline qpragma::hhl::observables::Observable<SIZE>::Observable(
    const std::vector<PauliTerm<SIZE>> & pauli_terms
): m_pauli_terms(pauli_terms) {}


template <uint64_t SIZE>
inline uint64_t qpragma::hhl::observables::Observable<SIZE>::size() const {
    return m_pauli_terms.size();
}


template <uint64_t SIZE>
inline std::vector<qpragma::hhl::observables::PauliTerm<SIZE>>::const_iterator
qpragma::hhl::observables::Observable<SIZE>::begin() const {
    return m_pauli_terms.cbegin();
}


template <uint64_t SIZE>
inline std::vector<qpragma::hhl::observables::PauliTerm<SIZE>>::const_iterator
qpragma::hhl::observables::Observable<SIZE>::end() const {
    return m_pauli_terms.cend();
}


template <uint64_t SIZE>
inline qpragma::hhl::observables::Observable<SIZE> &
qpragma::hhl::observables::Observable<SIZE>::operator*=(double constant) {
    (*this) = constant * (*this);
    return (*this);
}


// Operator
template <uint64_t SIZE>
inline qpragma::hhl::observables::Observable<SIZE> operator*(
    double constant, const qpragma::hhl::observables::Observable<SIZE> & obs
) {
    using namespace qpragma::hhl::observables;

    // Create new terms
    std::vector<PauliTerm<SIZE>> new_terms;
    new_terms.reserve(obs.size());

    // Update new terms
    std::transform(obs.begin(), obs.end(), std::back_inserter(new_terms),
                   [constant](auto term){ return PauliTerm<SIZE>(term.pauli_string(), constant * term.coeff()); });

    // Return result
    return Observable<SIZE>(new_terms);
}
