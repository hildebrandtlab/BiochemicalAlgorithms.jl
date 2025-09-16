# All on Handling molecules


## How can I create a simple molecule?

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

## How can I determine the element of an atom (C, N, …)?

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

## How can I identify backbone atoms?

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

## How can I get the one-letter (or three-letter) code out of a pdb file?

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
for chain in chains(sys)
    # get all residues from the current chain
    res = residues(chain)
    println("> Chain $(chain.name)")
    println("One-letter code:")
    println(join(one_letter_code.(res.name), ""))
    println("Three-letter code:")
    println(join(three_letter_code.(res), " "))
end
```

    > Chain E
    One-letter code:
    IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN
    Three-letter code:
    ILE VAL GLY GLY TYR THR CYS GLY ALA ASN THR VAL PRO TYR GLN VAL SER LEU ASN SER GLY TYR HIS PHE CYS GLY GLY SER LEU ILE ASN SER GLN TRP VAL VAL SER ALA ALA HIS CYS TYR LYS SER GLY ILE GLN VAL ARG LEU GLY GLU ASP ASN ILE ASN VAL VAL GLU GLY ASN GLU GLN PHE ILE SER ALA SER LYS SER ILE VAL HIS PRO SER TYR ASN SER ASN THR LEU ASN ASN ASP ILE MET LEU ILE LYS LEU LYS SER ALA ALA SER LEU ASN SER ARG VAL ALA SER ILE SER LEU PRO THR SER CYS ALA SER ALA GLY THR GLN CYS LEU ILE SER GLY TRP GLY ASN THR LYS SER SER GLY THR SER TYR PRO ASP VAL LEU LYS CYS LEU LYS ALA PRO ILE LEU SER ASP SER SER CYS LYS SER ALA TYR PRO GLY GLN ILE THR SER ASN MET PHE CYS ALA GLY TYR LEU GLU GLY GLY LYS ASP SER CYS GLN GLY ASP SER GLY GLY PRO VAL VAL CYS SER GLY LYS LEU GLN GLY ILE VAL SER TRP GLY SER GLY CYS ALA GLN LYS ASN LYS PRO GLY VAL TYR THR LYS VAL CYS ASN TYR VAL SER TRP ILE LYS GLN THR ILE ALA SER ASN
    > Chain I
    One-letter code:
    RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA
    Three-letter code:
    ARG PRO ASP PHE CYS LEU GLU PRO PRO TYR THR GLY PRO CYS LYS ALA ARG ILE ILE ARG TYR PHE TYR ASN ALA LYS ALA GLY LEU CYS GLN THR PHE VAL TYR GLY GLY CYS ARG ALA LYS ARG ASN ASN PHE LYS SER ALA GLU ASP CYS MET ARG THR CYS GLY GLY ALA
    > Chain E
    One-letter code:

    Three-letter code:

    > Chain I
    One-letter code:

    Three-letter code:

## How can I pick one single chain out of a system containing several chains?

This is often needed when a receptor and a ligand are co-complexed and you want to treat them separately.

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))

all_chains = chains(sys)
```

| **\#** | **idx** | **name** |
|-------:|:--------|:---------|
|      1 | 2       | E        |
|      2 | 1855    | I        |
|      3 | 2368    | E        |
|      4 | 2615    | I        |

This snippet will create separate PDB files for the two chains of the system:
Note: Currently fails see [issue 209](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/issues/209)

``` julia
#write_pdb("2ptc_chainE.pdb", all_chains[1])
#write_pdb("2ptc_chainI.pdb", all_chains[2])
```

## How can I map two configurations of the same protein onto each other?

``` julia
# read in the first protein
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
mol = first(molecules(sys))

# we will generate a second version of the protein by moving it around
sys2 = deepcopy(sys)
mol2 = first(molecules(sys2))
translate!(mol2, Vector3{Float32}(2.0, 1.0, 2.0))

# let's see how far the structures are apart
println("RMSD before mapping:\t", compute_rmsd(mol2, mol))

# now we have two proteins or system we can map together
map_rigid!(mol2, mol)

# let's see how far the structures are apart afterwards
println("RMSD after mapping:\t", compute_rmsd(mol2, mol))
```

    RMSD before mapping:    3.0
    RMSD after mapping: 2.2895498e-5

## How can I rotate an entire molecule?

``` julia
# read in the first protein
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
mol = first(molecules(sys))

v = Vector3{Float32}(0, 0, 0) # no translation
m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0) # counter clockwise rotation by 90 degree
r = RigidTransform(m, v)

# perform the transformation
rigid_transform!(mol, r)
```

    BiochemicalAlgorithms.Molecule{Float32}: (idx = 1, name = "COMPLEX (PROTEINASE/INHIBITOR)")

## How can I remove water molecules from a system?

``` julia
sys = load_pdb(ball_data_path("../test/data/1tgh.pdb"))

println("Number of atoms before removing water: ", natoms(sys))

# find all water fragments
ft = filter(frag -> frag.name == "HOH", fragments(sys))

# delete the found fragments, including all atoms and bonds
delete!(ft)

println("Number of atoms after removing water: ", natoms(sys))
```

    Number of atoms before removing water: 2355
    Number of atoms after removing water: 2301

# How can I identify atoms in a certain spatial proximity efficiently?

Identifying atoms in a certain spatial proximity can be done with `CellListMap.jl` which is already a dependency of BiochemicalAlgorithms.jl. Here is an example with a cutoff distance of 1.5:

``` julia
using CellListMap

sys = load_pdb(ball_data_path("../test/data/1tgh.pdb"))

neighbors = neighborlist([a.r for a in atoms(sys)], 1.5)
```

    1693-element Vector{Tuple{Int64, Int64, Float64}}:
     (1, 19, 0.9589779473125072)
     (2, 1, 1.4388685066986233)
     (3, 4, 1.4438002160243146)
     (4, 8, 1.4110367235570076)
     (5, 6, 1.423614076686212)
     (8, 9, 1.4528402026870852)
     (9, 10, 1.3692218168573846)
     (10, 11, 1.2292563742786062)
     (10, 12, 1.3467147068968133)
     (13, 12, 1.328092566103565)
     ⋮
     (2344, 2345, 0.96483334318895)
     (2344, 2346, 0.9578454673076238)
     (2347, 2348, 0.9646835173520671)
     (2347, 2349, 0.9621348854276189)
     (2350, 2351, 0.9633766194458323)
     (2352, 2350, 0.9584511836064512)
     (2353, 2354, 0.9596359268143725)
     (2353, 2355, 0.9774217099578173)
     (2354, 2355, 1.4942598847835975)

`neighborlist` returns a tuple consisting of indices for the neighbors and the computed distance between them. For more advanced use cases have a look at the [`CellListMap.jl` documentation](https://m3g.github.io/CellListMap.jl/stable/neighborlists/) .
