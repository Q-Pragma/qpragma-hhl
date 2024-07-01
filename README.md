# C++ Q-Pragma hybrid HHL
This repository proposes an implementation of **hybrid HHL** using the C++ hybrid quantum framework **Q-Pragma**. This framework is meant to be used in **HPC environments** and is inspired by **OpenMP**. It illustrates some features available in both OpenMP and Q-Pragma crucial for hybrid computing but also some Q-Pragma specific ones that would make OpenMP compatible with quantum-hybrid computing. This implementation uses algorithms from [1] and [2].

## Hybrid-HHL algorithm
The hybrid HHL algorithm is an altered version of the original HHL algorithm [4]. It uses a **hybrid pre-processing** of the data to provide a simpler quantum circuit solving the linear system of equations. The algorithm can be described as follows:
1. **Sampling** the eigenvalues of the problem using the Quantum Phase Estimation (QPE).
2. **Classical processing** of the eigenvalues to find similarities between their binary encoding.
3. Generate the HHL quantum circuit **exploiting** the previous processing to **reduce** the number of quantum operations needed.

## Q-Pragma implementation of hybrid HHL
The file hybrid_hhl.cpp provides an implementation of hybrid HHL. 

### Instruction offloading
This implementation mixes classical, hybrid quantum-classical, and pure quantum code. This multi-leveling requires two types of **offloading**:

- Offloading of **hybrid quantum-classical code** to the QPU controller: this is done using the *pragma quantum scope* directive. The code can be offloaded with both quantum and classical instructions.
- Offloading of **pure quantum code** to the QPU: this is done using the *pragma quantum routine*. The code offloaded is pure quantum, so only pure quantum instructions are usable.

### Classical control flow on quantum instructions
We provide an example of **classical control flow** on quantum instructions (see do-while loop, lines 236-246). This pattern can be useful in hybrid quantum-classical computing, for example when we want to select on a specific quantum state. **Post-selection** is a hybrid operation, relying on repeated measurements until the desired outcome is obtained.

### Quantum computing specific pragma directives
Q-Pragma also provides two directives to help develop quantum code. These rely on the **controllability** and **reversibility** of quantum operations.
- *pragma quantum ctrl*: **quantum controls** the following quantum instructions.
- *pragma quantum compute*: **compute** the following quantum instructions and then, when the end of the current scope is reached, automatically **uncompute** (apply the reverse) these insturctions.

## Conclusion
This implementation shows how to integrate quantum computing in current HPC programming environments. In [3], we compare Q-Pragma and OpenMP features and suggest ways to combine them in OpenMP.

## References
[1] Zhang, M., Dong, L., Zeng, Y. et al. Improved circuit implementation of the HHL algorithm and its simulations on QISKIT. Sci Rep 12, 13287 (2022). https://doi.org/10.1038/s41598-022-17660-8

[2] Lee, Y., Joo, J. & Lee, S. Hybrid quantum linear equation algorithm and its experimental test on IBM Quantum Experience. Sci Rep 9, 4778 (2019). https://doi.org/10.1038/s41598-019-41324-9

[3] K., G., R., A., R., C. B, Pragma-Based Approach for Hybrid Quantum-HPC Programming, in submission, (2024)

[4] Aram W. Harrow, Avinatan Hassidim, Seth Lloyd. Quantum algorithm for solving linear systems of equations, Physical Review Letters, (2009). http://dx.doi.org/10.1103/PhysRevLett.103.150502
