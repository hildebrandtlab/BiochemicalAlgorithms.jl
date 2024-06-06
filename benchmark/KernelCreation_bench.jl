
kernel_creation = SUITE["Kernel"]["Creation"] 

n = 40000
sys = System()
kernel_creation["AtomCreation"] = @benchmarkable [Atom($sys,i, Elements.C) for i in 1:$n]


sys1 = System()
mol = Molecule(sys1)
chain = Chain(mol)
kernel_creation["FragmentCreation"] = @benchmarkable [Fragment($chain,i) for i in 1:$n]


sys2 = System()
mol = Molecule(sys2)
chain2 = Chain(mol)
kernel_creation["ResidueCreation"] = @benchmarkable [Residue($chain, i, AminoAcid('D')) for i in 1:$n]

sys3 = System()
kernel_creation["MoleculeCreation"] = @benchmarkable [Molecule($sys3) for i in 1:$n]

kernel_creation["SystemCreation"] = @benchmarkable [System() for i in 1:$n]
