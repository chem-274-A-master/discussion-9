# Input provenance

`check.xyz` contains four synthetic points, not a chemically realistic structure. Its unique closest pair is indices 1 and 3, separated by 2 angstrom.

`caffeine.sdf` and `caffeine.xyz` were generated for this exercise from caffeine SMILES `Cn1c(=O)c2c(ncn2C)n(C)c1=O`, using RDKit embedding with seed 274 followed by MMFF optimization. They are teaching fixtures, not downloaded PubChem records or experimental coordinates. RDKit is not required to use them. Both contain 24 atoms, including 10 hydrogens.

`../python_debugging/caffeine_reference.png` was generated from the SDF with the corrected plotting implementation. It shows the 14 non-hydrogen atoms and 15 bonds. A graph layout is not a depiction of the 3D coordinates; rotations, reflection and different placement are acceptable.
