# Getting started


BiochemicalAlgorithms.jl\`s design is inspired by its predecessor [`BALL`](https://github.com/ball-project/ball). BALL and BiochemicalAlgorithms.jl both have a set of KERNEL classes or core classes to represent molecular entities such as Molecules, Proteins, Atoms, Bonds etc. The following class diagramm shows the different classes and their relationship to each others.

![UML class diagram *Note*: Only the most important functionalities are shown here.](uml.png)

UML class diagram *Note*: Only the most important functionalities are shown here.

In the center of all classes resides the [`System`](@ref). So let’s see it in action with a simple peptide:

``` julia
using BiochemicalAlgorithms

s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
s
```

    System{Float32}: AlaAla.pdb
      23 atoms
       0 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments

This will give us an overview of the objects associated with our system. It contains 23 atoms, 0 bonds, 1 molecule, 1 chain, 0 secondary structures and 2 fragments. We look at each group separately in the following sections.

If we do not explicitely create a system, the data will be stored in the default system:

``` julia
load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
```

    System{Float32}: AlaAla.pdb
      23 atoms
       0 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments

## Atom

We can access the atoms like this:

``` julia
sys_atms = atoms(s)

natoms(sys_atms)
```

    23

Let’s play around with it:

``` julia
println("Atom elements:")
for a in sys_atms
    print(a.element)
end
println()

n1= atom_by_name(s, "N")

println(n1.idx)
```

    Atom elements:
    NCCOHHHHHHCHNCCHHHHCHOO
    5

We can access an atom by its name, which will return the first atom matching the given name. If there are more atoms with the same name and we are not interested in the first one, we can access the atom by the atom index.
The atom index is unique inside a system but not necessarily starting by one.

``` julia
println("Atom ids:")
for a in sys_atms
    println(a.idx, " ", a.element)
end

n2 = atom_by_idx(s, 17)
```

    Atom ids:
    5 N
    6 C
    7 C
    8 O
    9 H
    10 H
    11 H
    12 H
    13 H
    14 H
    15 C
    16 H
    17 N
    18 C
    19 C
    20 H
    21 H
    22 H
    23 H
    24 C
    25 H
    26 O
    27 O

    BiochemicalAlgorithms.Atom{Float32}: (idx = 17, number = 13, element = BiochemicalAlgorithms.Elements.N, name = "N", atom_type = "", r = Float32[0.133, 0.98, 2.268], v = Float32[0.0, 0.0, 0.0], F = Float32[0.0, 0.0, 0.0], formal_charge = 0, charge = 0.0f0, radius = 0.0f0)

[`Atom`](@ref) features the atom number as well. The atom number does NOT have to be unique inside the system!

``` julia
for a in sys_atms
    println(a.idx, " ", a.number)
end
```

    5 1
    6 2
    7 3
    8 4
    9 5
    10 6
    11 7
    12 8
    13 9
    14 10
    15 11
    16 12
    17 13
    18 14
    19 15
    20 16
    21 17
    22 18
    23 19
    24 20
    25 21
    26 22
    27 23

Next, we will have a look at `Bond`.

## Bond

BiochemicalAlgorithms.jl supports different BondTypes:

``` julia
BiochemicalAlgorithms.BondOrderType
```

    Enum type BondOrder.T <: Enum{Int32} with 6 instances:
     BondOrder.Single    = 1
     BondOrder.Double    = 2
     BondOrder.Triple    = 3
     BondOrder.Quadruple = 4
     BondOrder.Aromatic  = 50
     BondOrder.Unknown   = 100

Typically, bonds are not included in the PDB-FileFormat and have to be computed. This can be achieved with the FragmentDataBase, which contains known fragments (all amino acids such as alanine, all nucleotides, same ions …). If you are interested in the fragments contained in the default FragmentDB, take a look at data/fragments of the repository. There you find all the fragments, each one stored in a separate json format.
So let’s take a look at the FragmentDB and bonds:

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
nbonds(s)
```

    0

So far, the system has no bonds. We can also visualize the structure using BiochemicalVisualization.jl.

If BiochemicalVisualization is not yet installed, you can add it via:

``` julia
using Pkg
Pkg.add("BiochemicalVisualization")
```

       Resolving package versions...
          Compat entries added for BiochemicalVisualization
        Updating `~/local/BiochemicalAlgorithms.jl/Project.toml`
      [00dfb448] + BiochemicalVisualization v0.3.2
      No Changes to `~/local/BiochemicalAlgorithms.jl/Manifest.toml`
    ┌ Warning: Circular dependency detected.
    │ Precompilation will be skipped for dependencies in this cycle:
    │  ┌ BiochemicalAlgorithms
    │  └─ BiochemicalVisualization
    └ @ Base.Precompilation precompilation.jl:650

``` julia
using BiochemicalVisualization

# ball_and_stick(s)  uncomment for visualization
```

Now, let’s create some bonds:

``` julia
fdb = FragmentDB() # default FragmentDB
normalize_names!(s, fdb) # in case our input PDB file uses a strange naming standard
reconstruct_fragments!(s, fdb) # in case our input file misses some atoms
build_bonds!(s, fdb) # create the bonds
# ball_and_stick(s) uncomment for visualization
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 22 bonds

And finally, we do have our bonds. This example demonstrate the importance of visualizing for the understanding of structural data.
Note: The exclamation mark behind the function name indicates that the system will be changed through the functions. The system now contains bonds and the residues carries flags:

``` julia
println(s)
for i in residues(s)
    println(i.flags)
end
```

    BiochemicalAlgorithms.System{Float32} with 23 atoms and 22 bonds (AlaAla.pdb)
    Set([:N_TERMINAL])
    Set([:C_TERMINAL])

Let’s take a deeper look at the bonds:

``` julia
for bond in bonds(s)
    println(bond)
end
```

    BiochemicalAlgorithms.Bond{Float32}: (idx = 28, a1 = 5, a2 = 13, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 29, a1 = 5, a2 = 6, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 30, a1 = 5, a2 = 9, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 31, a1 = 5, a2 = 11, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 32, a1 = 6, a2 = 16, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 33, a1 = 6, a2 = 7, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 34, a1 = 6, a2 = 15, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 35, a1 = 7, a2 = 8, order = BiochemicalAlgorithms.BondOrder.Double)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 36, a1 = 10, a2 = 15, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 37, a1 = 12, a2 = 15, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 38, a1 = 14, a2 = 15, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 39, a1 = 17, a2 = 20, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 40, a1 = 17, a2 = 18, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 41, a1 = 18, a2 = 25, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 42, a1 = 18, a2 = 19, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 43, a1 = 18, a2 = 24, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 44, a1 = 19, a2 = 26, order = BiochemicalAlgorithms.BondOrder.Double)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 45, a1 = 19, a2 = 27, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 46, a1 = 21, a2 = 24, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 47, a1 = 22, a2 = 24, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 48, a1 = 23, a2 = 24, order = BiochemicalAlgorithms.BondOrder.Single)
    BiochemicalAlgorithms.Bond{Float32}: (idx = 49, a1 = 7, a2 = 17, order = BiochemicalAlgorithms.BondOrder.Single)

Bond are only formed inside the same system and not across systems. The atoms forming the bonds can be accessed through the following:

``` julia
for bond in bonds(s)
    println("Bond between atom ", bond.a1, " and ", bond.a2, " and has an order: ", bond.order)
end
```

    Bond between atom 5 and 13 and has an order: Single
    Bond between atom 5 and 6 and has an order: Single
    Bond between atom 5 and 9 and has an order: Single
    Bond between atom 5 and 11 and has an order: Single
    Bond between atom 6 and 16 and has an order: Single
    Bond between atom 6 and 7 and has an order: Single
    Bond between atom 6 and 15 and has an order: Single
    Bond between atom 7 and 8 and has an order: Double
    Bond between atom 10 and 15 and has an order: Single
    Bond between atom 12 and 15 and has an order: Single
    Bond between atom 14 and 15 and has an order: Single
    Bond between atom 17 and 20 and has an order: Single
    Bond between atom 17 and 18 and has an order: Single
    Bond between atom 18 and 25 and has an order: Single
    Bond between atom 18 and 19 and has an order: Single
    Bond between atom 18 and 24 and has an order: Single
    Bond between atom 19 and 26 and has an order: Double
    Bond between atom 19 and 27 and has an order: Single
    Bond between atom 21 and 24 and has an order: Single
    Bond between atom 22 and 24 and has an order: Single
    Bond between atom 23 and 24 and has an order: Single
    Bond between atom 7 and 17 and has an order: Single

Bonds can be accessed via their respective index, deleted, and be put in the system

``` julia
bond = bond_by_idx(s, 28)
println("Bond index: ", bond.idx, ", between atom ", bond.a1, " and ", bond.a2, ", with order: ", bond.order)

println("Bonds in the system: ", nbonds(s))

delete!(bond)

println("Bonds in the system after deletion: ", nbonds(s))

bond = Bond(s, 5,13, BondOrder.Single)

println("Bonds in the system after deletion: ", nbonds(s))
```

    Bond index: 28, between atom 5 and 13, with order: Single
    Bonds in the system: 22
    Bonds in the system after deletion: 21
    Bonds in the system after deletion: 22

Be careful: The bond index is now different because it is a new Bond object:

``` julia
println("Bond index: ", bond.idx, ", between atom ", bond.a1, " and ", bond.a2, ", with order: ", bond.order)
```

    Bond index: 50, between atom 5 and 13, with order: Single

Another helpful functionalities when it comes to bonds are the following:

``` julia
# get the atoms sharing a bond
a1, a2 = get_partners(bond)

# or if we have only given one atom
atom = atom_by_idx(s,12)
bond = bonds(atom)[1]
partner_atom = get_partner(bond, atom)
println(partner_atom.idx)
println("The bond length is ", bond_length(bond))
```

    15
    The bond length is 1.0907414

Sometimes it is necessary to discriminate between hydrogen and nonhydrogen_bonds. Since this a common use case, there are built in functionalities:

## Fragment

Fragments represent molecule fragments. In context of proteins or nucleic acids, the fragment can be used to represent residues or nucleotides. See the documentation for the definition of the [`FragmentVariantType`](@ref)

``` julia
nfragments(s)

for frag in fragments(s)
    println(frag.name, " ", frag.idx)
end
n2 = atom_by_idx(s, 5)

println("The N atom belongs to residue: ", n2.fragment_idx)
```

    ALA 3
    ALA 4
    The N atom belongs to residue: 3

We can get the fragment as well directly:

``` julia
fragment = parent_fragment(n2)
println(fragment)
res2 = fragment_by_idx(s, 3)
```

    BiochemicalAlgorithms.Fragment{Float32}: (idx = 3, number = 1, name = "ALA")

    BiochemicalAlgorithms.Fragment{Float32}: (idx = 3, number = 1, name = "ALA")

Just like atom indices, fragment indices not necessarily follow a certain order but are unique inside the same system object.

``` julia
res1 = fragments(s)[1]
res1.variant # 1 = unknown, 2 = residue, 3 = nucleotide
```

    FragmentVariant.Residue = 2

Alternatively, it is also possible to check if this is a nucleotide:

``` julia
println(is_nucleotide(res1))

println("Number of fragments before push: ", nfragments(s))
push!(chain_by_idx(s, 2), res1)
println("Number of fragments after push: ", nfragments(s))
delete!(res1)
println("Number of fragments after deletion: ", nfragments(s))

for frag in fragments(s)
    println(frag.name, " ", frag.idx)
end
s
```

    false
    Number of fragments before push: 2
    Number of fragments after push: 3
    Number of fragments after deletion: 2
    ALA 4
    ALA 51

    System{Float32}: AlaAla.pdb
      11 atoms
      10 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments

Or a residue:

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb")) # restore original data
res1 = fragments(s)[1]
println("is residue ", isresidue(res1))
```

    is residue true

There are more functionalities related to fragments. Have a look at the documentation for [`Fragment`](@ref) for more details.

## Residue

As described in the Fragments section, [`Residue`](@ref) is a Fragment with `FragmentVariant.Residue` and typically describes an amino acid:

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb")) # restore original data
res = fragments(s)[1]
println("is_amino_acid(): ", is_amino_acid(res))
println("is_c_terminal(): ", is_c_terminal(res))
println("is_n_terminal(): ", is_n_terminal(res))
```

    is_amino_acid(): true
    is_c_terminal(): false
    is_n_terminal(): false

Some functionality is only available after preprocessing through the FragmentDB:

``` julia
fdb = FragmentDB()
normalize_names!(s, fdb)
reconstruct_fragments!(s, fdb)
build_bonds!(s, fdb)
println("is_c_terminal(): ", is_c_terminal(res))
println("is_n_terminal(): ", is_n_terminal(res))
res2 = fragments(s)[2]
println("is_c_terminal(): ", is_c_terminal(res2))
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 22 bonds
    is_c_terminal(): false
    is_n_terminal(): true
    is_c_terminal(): true

## Nucleotide

[`Nucleotide`](@ref) is a Fragment with `FragmentVariant.Nucleotide`:

``` julia
s = System()
chain = Chain(Molecule(s))
n1 = Nucleotide(chain, 1;
            name = "my nucleotide",
            properties = Properties(:first => 'a', :second => "b"),
            flags = Flags([:third])
        )
        n2 = Nucleotide(chain, 2)
```

    BiochemicalAlgorithms.Fragment{Float32}: (idx = 4, number = 2, name = "")

## Molecule

Molecules are used to distinguish between proteins and non-proteins (see [`Molecules`](@ref)

``` julia
s = load_pdb(ball_data_path("../test/data/5PTI.pdb"))
m = molecules(s)
m
```

| **\#** | **idx** | **name** |
|-------:|:--------|:---------|
|  **1** | 1       | 5PTI.pdb |

Like in the other cases before, we get a table of molecules. Which we can access similarly to atoms and bonds:

``` julia
println(typeof(m))

molecule_by_idx(s, 1)
println("The number of molecules in the system: ", nmolecules(s))
println("The number of proteins in the system: ", nproteins(s))
```

    MoleculeTable{Float32}
    The number of molecules in the system: 1
    The number of proteins in the system: 0

Molecules have a variant which is `MoleculareVariant.None` by default. Users can decide, if a molecule is to be considered as a protein. `isprotein` checks if the `MoleculeVariant.Protein` is set:

``` julia
println("is a protein: ", isprotein(m[1]))
fdb = FragmentDB()
normalize_names!(s, fdb)
reconstruct_fragments!(s, fdb)
build_bonds!(s, fdb)

println("is a protein: ", isprotein(m[1]))

m[1].variant = MoleculeVariant.Protein

println("is a protein: ", isprotein(m[1]))
```

    is a protein: false
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for PO4:70
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:80
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:101
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:102
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:105
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:110
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:111
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:112
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:113
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:116
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:117
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:119
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:121
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:122
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:125
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:126
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:127
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:129
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:133
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:134
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:138
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:140
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:143
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:144
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:145
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:146
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:156
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:157
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:158
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:159
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:160
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:200
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:201
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:202
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:203
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:204
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:205
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:209
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:210
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:211
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:212
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:214
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:216
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:217
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:218
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:219
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:220
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:223
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:225
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:302
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:304
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:310
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:311
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:312
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:313
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:314
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:315
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:316
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:317
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:318
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:319
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:320
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:321
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for DOD:322
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    ┌ Warning: reconstruct_fragments!(): could not find reference fragment for UNX:324
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/reconstruct_fragments.jl:177
    [ Info: reconstruct_fragments!(): added 109 atoms.
    ┌ Warning: build_bonds!(): could not find reference fragment for PO4.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for DOD.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    ┌ Warning: build_bonds!(): could not find reference fragment for UNX.
    └ @ BiochemicalAlgorithms ~/local/BiochemicalAlgorithms.jl/src/preprocessing/build_bonds.jl:14
    [ Info: build_bonds!(): built 912 bonds
    is a protein: false
    is a protein: true

Similar to molecules we can access proteins:

``` julia
println("Number of proteins in the system: ", proteins(s))
pti =protein_by_idx(s,1)
```

    Number of proteins in the system: BiochemicalAlgorithms.MoleculeTable{Float32} with 1 rows

    BiochemicalAlgorithms.Molecule{Float32}: (idx = 1, name = "5PTI.pdb")

A very useful functionality is the access the parent molecule or protein a specific atom is belonging to in cases where we have several molecules to deal with.

``` julia
atom = atom_by_idx(s,126)
println("This atom belongs to protein: ", parent_protein(atom))
```

    This atom belongs to protein: BiochemicalAlgorithms.Molecule{Float32}: (idx = 1, name = "5PTI.pdb")

## Chains

Chains can be considered as an aggregation of fragments (either residues or nucleotides).

``` julia

s = load_pdb(ball_data_path("../test/data/4hhb.pdb"))

println("Number of chains in the system: ", nchains(s))

for chain in chains(s)
    println("Chain index: ", chain.idx)
end

chain = chain_by_idx(s, 2)
println("Chain 2 contains ", nfragments(chain), " fragments.")

println("The parent molecule of chain 2 is: ", parent_molecule(chain))
```

    Number of chains in the system: 4
    Chain index: 2
    Chain index: 201
    Chain index: 407
    Chain index: 609
    Chain 2 contains 198 fragments.
    The parent molecule of chain 2 is: BiochemicalAlgorithms.Molecule{Float32}: (idx = 1, name = "4hhb.pdb")

## Congratulation

You are now familiar with the most important core entities of BiochemicalAlgorithms.jl
