# All on Handling molecules

``` julia
using BiochemicalAlgorithms
```

# How can I create a simple molecule?

``` julia
# create a system first
sys = System{Float32}() # this system will be of single precision, i.e., atom positions, velocities... 

# as well as a molecule
h2o = Molecule(sys)

# create system atoms
o1 = Atom(h2o, 1, Elements.O)
h1 = Atom(h2o, 2, Elements.H)
h2 = Atom(h2o, 3, Elements.H)

# set positions of the atoms
# o1.r = [0, 0, 0]  <-- this is the default value!
h1.r = [1, 0, 0]
h2.r = [cos(deg2rad(105)), sin(deg2rad(105)), 0)

# add bonds
Bond(h2o, o1.idx, h1.idx, BondOrder.Single)
Bond(h2o, o1.idx, h2.idx, BondOrder.Single)

println("Number of atoms: ", natoms(h2o))
println("Number of bonds: ", nbonds(h2o))
```

# How can I determine the element type of an atom (i.e., C,N,…)?

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

for atom in atoms(s)
    if atom.element == Elements.C
        println("This is a carbon atom!")
    end
end
```

# How can I identify backbone atoms?

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

for atom in atoms(s)
    println(atom.element)
    if atom.name in ["C", "O", "N", "HA"]
        println("This is a backbone atom!")
    end
end
```

# How can I get the one-letter code out of a PDB file?

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
for chain in chains(sys)
    print.(res.type for res in residues(sys))
end
```

# How can I pick one single chain out of a system containing several chains?

This is often needed when a receptor and a ligand are co-complexed and
you want to treat them separately. This snippet will create copies of
the first two chains and strores them in separate pdb files.

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
chainA = deepcopy(chains(sys)[1])
chainB = deepcopy(chains(sys)[2])

write_pdb("2ptc_chainA.pdb", chainA)
write_pdb("2ptc_chainB.pdb", chainB)
```

# How can I map two configurations of the same protein onto each other?

``` julia
# read in the first protein
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
mol = first(molecules(sys))

# we will generate a second version of the protein by moving it around
sys2 = deepcopy(sys)
mol2 = first(molecules(sys2))
translate!(mol2, Vector3{Float32}(2.0,1.0,2.0))

# let's see how far the structures are apart
println(compute_rmsd(mol2, mol))

# now we have two proteins or system we can map together
map_rigid!(mol2, mol)

# let's see how far the structures are apart afterwards
println(compute_rmsd(mol2, mol))
```
