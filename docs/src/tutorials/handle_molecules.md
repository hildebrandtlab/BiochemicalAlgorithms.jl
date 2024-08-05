# All on Handling molecules


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
h2.r = [cos(deg2rad(105)), sin(deg2rad(105)), 0]

# add bonds
Bond(h2o, o1.idx, h1.idx, BondOrder.Single)
Bond(h2o, o1.idx, h2.idx, BondOrder.Single)

println("Number of atoms: ", natoms(h2o))
println("Number of bonds: ", nbonds(h2o))
```

    Number of atoms: 3
    Number of bonds: 2

# How can I determine the element of an atom (C, N, …)?

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
for atom in atoms(sys)
    println("Atom no.: $(atom.number), element: $(atom.element)")
end
```

    Atom no.: 1, element: N
    Atom no.: 2, element: C
    Atom no.: 3, element: C
    Atom no.: 4, element: O
    Atom no.: 5, element: H
    Atom no.: 6, element: H
    Atom no.: 7, element: H
    Atom no.: 8, element: H
    Atom no.: 9, element: H
    Atom no.: 10, element: H
    Atom no.: 11, element: C
    Atom no.: 12, element: H
    Atom no.: 13, element: N
    Atom no.: 14, element: C
    Atom no.: 15, element: C
    Atom no.: 16, element: H
    Atom no.: 17, element: H
    Atom no.: 18, element: H
    Atom no.: 19, element: H
    Atom no.: 20, element: C
    Atom no.: 21, element: H
    Atom no.: 22, element: O
    Atom no.: 23, element: O

You can also filter for specific elements:

``` julia
c_atoms = filter(atom -> atom.element == Elements.C, atoms(sys))
for atom in c_atoms
    println("Atom no.: $(atom.number), element: $(atom.element)")
end
```

    Atom no.: 2, element: C
    Atom no.: 3, element: C
    Atom no.: 11, element: C
    Atom no.: 14, element: C
    Atom no.: 15, element: C
    Atom no.: 20, element: C

# How can I identify backbone atoms?

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
for atom in atoms(sys)
    print(atom.element)
    if atom.name in ["C", "O", "N", "HA"]
        print(" <-- This is a backbone atom!")
    end
    println()
end
```

    N <-- This is a backbone atom!
    C
    C <-- This is a backbone atom!
    O <-- This is a backbone atom!
    H
    H
    H
    H
    H
    H
    C
    H <-- This is a backbone atom!
    N <-- This is a backbone atom!
    C
    C <-- This is a backbone atom!
    H
    H
    H
    H
    C
    H <-- This is a backbone atom!
    O <-- This is a backbone atom!
    O

# How can I get the one-letter code out of a pdb file?

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
for chain in chains(sys)
    # get all residues from the current chain
    res = residues(chain).name
    println("> Chain $(chain.name)")
    println(join(one_letter_code.(res), ""))
end
```

    > Chain E
    IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAYGLEGKGDSCQGDSGGPVVCSGKLQGIVSWGSGCQAKNKPGVYTKVCNYVSWIKQTIASN
    > Chain I
    RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA

# How can I pick one single chain out of a system containing several chains?

This is often needed when a receptor and a ligand are co-complexed and you want to treat them separately.

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))

all_chains = chains(sys)
```

|  \# | idx | name |
|----:|:----|:-----|
|   1 | 2   | E    |
|   2 | 349 | I    |

BiochemicalAlgorithms.ChainTable{Float32} with 2 rows:

This snippet will create separate PDB files for the two chains of the system:

``` julia
write_pdb("2ptc_chainE.pdb", all_chains[1])
write_pdb("2ptc_chainI.pdb", all_chains[2])
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
println("RMSD before mapping:\t", compute_rmsd(mol2, mol))

# now we have two proteins or system we can map together
map_rigid!(mol2, mol)

# let's see how far the structures are apart afterwards
println("RMSD after mapping:\t", compute_rmsd(mol2, mol))
```

    RMSD before mapping:    3.0
    RMSD after mapping: 3.881571e-6
