# OrbitalProjection
Projects one set of orbitals onto another for ChronusQuantum (CQ)

# Constraints
1. Same basis set and orbital saved to /SCF/MO1 in each bin file.
2. Assumes 2-component size for both MOs. For RSCF, ROSCF, and USCF, read your MOs into CQ as a scratch file with a 2-component reference and do SCF.GUESS=READMOS and SCF.ALG=SKIP.

# Usage
`python moProj.py binfile1.bin binfile2.bin`

The above assumes that binfile1.bin contains the MOs you will be using as a basis.
