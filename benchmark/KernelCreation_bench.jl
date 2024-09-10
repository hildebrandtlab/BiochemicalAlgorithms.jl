
kernel_creation = SUITE["Kernel"]["Creation"] 

sys = System()
kernel_creation["AtomCreation"] = @benchmarkable Atom($sys, 1, Elements.C) samples = 1000000


sys1 = System()
mol = Molecule(sys1)
chain = Chain(mol)
kernel_creation["FragmentCreation"] = @benchmarkable Fragment($chain,1) samples = 1000000


sys2 = System()
mol = Molecule(sys2)
chain2 = Chain(mol)

kernel_creation["ResidueCreation"] = @benchmarkable Residue($chain, 1, name = "ALA") samples = 1000000


sys3 = System()
kernel_creation["MoleculeCreation"] = @benchmarkable Molecule($sys3) samples = 1000000


kernel_creation["SystemCreation"] = @benchmarkable System() samples = 1000000
