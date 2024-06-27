
kernel_creation = SUITE["Kernel"]["Creation"] 

sys = System()
kernel_creation["AtomCreation"] = @benchmarkable Atom($sys,i, Elements.C)


sys1 = System()
mol = Molecule(sys1)
chain = Chain(mol)
kernel_creation["FragmentCreation"] = @benchmarkable Fragment($chain,i)


sys2 = System()
mol = Molecule(sys2)
chain2 = Chain(mol)
kernel_creation["ResidueCreation"] = @benchmarkable Residue($chain, i, AminoAcid('D'))


sys3 = System()
kernel_creation["MoleculeCreation"] = @benchmarkable Molecule($sys3)


kernel_creation["SystemCreation"] = @benchmarkable System()
