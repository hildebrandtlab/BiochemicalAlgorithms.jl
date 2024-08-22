# All on iteration


When working with molecular entities, we want to run over all atoms of a system, over all chains, etc. In this tutorial we will learn how this can be done.

## Molecular systems

In `BiochemicalAlgorithms.jl` atoms and bonds are existing inside a `System`. Typically, molecular data is stored in molecular data formats such as PDB. The latter can be directly read into a system.

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
```

    System with 23 atoms (AlaAla.pdb)

You can, e.g., print all atoms of the given system as a table:

``` julia
atoms(sys)
```

| \# | idx | number | element | name | atom_type | r | v | F | formal_charge | charge | radius |
|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 1 | 5 | 1 | N | N |  | Float32\[-1.45, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 2 | 6 | 2 | C | CA |  | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 3 | 7 | 3 | C | C |  | Float32\[0.495, 0.0, 1.437\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 4 | 8 | 4 | O | O |  | Float32\[1.235, -0.911, 1.838\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 5 | 9 | 5 | H | 1H |  | Float32\[-1.788, 0.918, 0.25\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 6 | 10 | 6 | H | 1HB |  | Float32\[1.642, 1.124, -0.864\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 7 | 11 | 7 | H | 2H |  | Float32\[-1.801, -0.26, -0.911\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 8 | 12 | 8 | H | 2HB |  | Float32\[0.154, 1.136, -1.827\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 9 | 13 | 9 | H | 3H |  | Float32\[-1.749, -0.665, 0.704\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 10 | 14 | 10 | H | 3HB |  | Float32\[0.258, 2.118, -0.351\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 11 | 15 | 11 | C | CB |  | Float32\[0.554, 1.175, -0.813\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 12 | 16 | 12 | H | HA |  | Float32\[0.341, -0.928, -0.46\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 13 | 17 | 13 | N | N |  | Float32\[0.133, 0.98, 2.268\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 14 | 18 | 14 | C | CA |  | Float32\[0.605, 0.98, 3.639\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 15 | 19 | 15 | C | C |  | Float32\[0.414, -0.408, 4.228\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 16 | 20 | 16 | H | H |  | Float32\[-0.476, 1.73, 1.938\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 17 | 21 | 17 | H | 1HB |  | Float32\[0.325, 2.106, 5.472\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 18 | 22 | 18 | H | 2HB |  | Float32\[0.065, 3.039, 3.988\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 19 | 23 | 19 | H | 3HB |  | Float32\[-1.16, 1.868, 4.519\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 20 | 24 | 20 | C | CB |  | Float32\[-0.089, 2.07, 4.463\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 21 | 25 | 21 | H | HA |  | Float32\[1.676, 1.185, 3.63\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 22 | 26 | 22 | O | O |  | Float32\[0.531, -1.358, 3.421\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 23 | 27 | 23 | O | OXT |  | Float32\[0.462, -0.512, 5.473\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |

BiochemicalAlgorithms.AtomTable{Float32} with 23 rows:

Single columns can be directly accessed by their name:

``` julia
atoms(sys).name
```

    23-element SystemComponentTableCol{String}:
     "N"
     "CA"
     "C"
     "O"
     "1H"
     "1HB"
     "2H"
     "2HB"
     "3H"
     "3HB"
     ⋮
     "C"
     "H"
     "1HB"
     "2HB"
     "3HB"
     "CB"
     "HA"
     "O"
     "OXT"

You can also directly query the number of atoms via `natoms`:

``` julia
natoms(sys)
```

    23

Similar functions exist for other system components, including bonds, molecules, chains, and fragments.

## How can I iterate over all atoms of a system?

We can just as easily iterate over all atoms:

``` julia
for atom in atoms(sys)
    # do something with this atom
end
```

Or, for individual columns:

``` julia
for atom_name in atoms(sys).name
    # do something with this atom name
end
```

## How can I iterate over specific atoms?

In many scenarios, we only want to iterate over a subset of atoms fulfilling a specific criteria. For example, here we only want to get the positions of the $C_\alpha$ atoms:

``` julia
ca_atoms = filter(atom -> atom.name == "CA", atoms(sys))
```

| \# | idx | number | element | name | atom_type | r | v | F | formal_charge | charge | radius |
|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 1 | 6 | 2 | C | CA |  | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 2 | 18 | 14 | C | CA |  | Float32\[0.605, 0.98, 3.639\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |

BiochemicalAlgorithms.AtomTable{Float32} with 2 rows:

And here we only want the heavy atoms:

``` julia
heavy_atoms = filter(atom -> atom.element != Elements.H, atoms(sys))
```

| \# | idx | number | element | name | atom_type | r | v | F | formal_charge | charge | radius |
|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 1 | 5 | 1 | N | N |  | Float32\[-1.45, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 2 | 6 | 2 | C | CA |  | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 3 | 7 | 3 | C | C |  | Float32\[0.495, 0.0, 1.437\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 4 | 8 | 4 | O | O |  | Float32\[1.235, -0.911, 1.838\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 5 | 15 | 11 | C | CB |  | Float32\[0.554, 1.175, -0.813\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 6 | 17 | 13 | N | N |  | Float32\[0.133, 0.98, 2.268\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 7 | 18 | 14 | C | CA |  | Float32\[0.605, 0.98, 3.639\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 8 | 19 | 15 | C | C |  | Float32\[0.414, -0.408, 4.228\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 9 | 24 | 20 | C | CB |  | Float32\[-0.089, 2.07, 4.463\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 10 | 26 | 22 | O | O |  | Float32\[0.531, -1.358, 3.421\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |
| 11 | 27 | 23 | O | OXT |  | Float32\[0.462, -0.512, 5.473\] | Float32\[0.0, 0.0, 0.0\] | Float32\[0.0, 0.0, 0.0\] | 0 | 0.0 | 0.0 |

BiochemicalAlgorithms.AtomTable{Float32} with 11 rows:

Instead of a system, we can use any other system component as an argument to the `atoms` function to only list the atoms of the same component. For example, the following system contains two chains:

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))
chains(sys)
```

|  \# | idx | name |
|----:|:----|:-----|
|   1 | 2   | E    |
|   2 | 349 | I    |

BiochemicalAlgorithms.ChainTable{Float32} with 2 rows:

In order to get all atoms of the first chain, we can use:

``` julia
chainE = first(chains(sys))
for atom in atoms(chainE)
    # do something with this atom of the first chain
end
```

We can even combine this approach with filtering to get, e.g., only the $C_\alpha$ atoms of the first chain:

``` julia
for atom in filter(atom -> atom.name == "CA", atoms(chainE))
    # do something with this CA atom of the first chain
end
```

## How can I iterate over all bonds of a system?

Bonds are not explicitely stored in the PDB format but are rather inferred after reading the data into a system using the fragment database `FragmentDB`:

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

# bonds are not contained in the pdb file
nbonds(sys)
```

    0

``` julia
# use the fragment database for normalizing naming schemas between molecular file formats, reconstruction of missing parts of the structure and building the bonds
fdb = FragmentDB()

normalize_names!(sys, fdb)
reconstruct_fragments!(sys, fdb)
build_bonds!(sys, fdb)

nbonds(sys)
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 22 bonds

    22

Similar to the atom iteration, we can iterate over all bonds of a sysem:

``` julia
for bond in bonds(sys)
    # do something with this bond
end
```

## How can I iterate over all fragments/residues/nucleotides of a system?

Chains usually contain fragments, e. g., residues or nucleotides and can either be queried collectively…

``` julia
sys = load_pdb(ball_data_path("../test/data/2ptc.pdb"))

nfragments(sys)
```

    439

…or by a given type:

``` julia
println("Number of residues: ", nresidues(sys))
println("Number of nucleotides: ", nnucleotides(sys))
```

    Number of residues: 281
    Number of nucleotides: 0

Just like in the examples above, the functions can also be used with other atom containers:

``` julia
nresidues.(chains(sys))
```

    2-element Vector{Int64}:
     223
      58
