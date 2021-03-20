# ComplexAdvancedHMC.jl

Extendes AdvancedHMC package to complex-valued parameters.

Changes from AdvancedHMC:
- General support for complex-valued parameters while still maintaining most of the support for real-valued parameters.
- Removed support for DenseMatrices, only supports uncorrelated noise.
- Added HermitianMetric to support hermitian parameters.

Dynamic sampling and stopping criterion not supported yet.
