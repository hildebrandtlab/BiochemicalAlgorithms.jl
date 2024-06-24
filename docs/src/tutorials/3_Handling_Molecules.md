# All on Handling molecules


``` julia
using BiochemicalAlgorithms
```

# How can I create a simple molecule?

``` julia
# create a system first

sys = System{Float32}() # this system will be of single precision, i.e., atom positions, velocities... 

h2o = Molecule(sys)
# create system atoms
o1 = Atom(h2o, 1, Elements.O)
h1 = Atom(h2o, 2, Elements.H)
h2 = Atom(h2o, 3, Elements.H)

o1.idx
# set positions of the atoms
# o1.r = Vector3{Float}(0, 0, 0)  <-- this is the default value!
h1.r = Vector3{Float32}(1, 0, 0)
h2.r = Vector3{Float32}(cos(105 * π / 180), sin(105 * π / 180), 0)

# add bonds
Bond(h2o, o1.idx, h1.idx, BondOrder.Single)
Bond(h2o, o1.idx, h2.idx, BondOrder.Single)
# 

println(natoms(h2o))
println(nbonds(h2o))
println.(atom for atom in atoms(h2o))
```

    3
    2
    Atom{Float32}: (idx = 2, number = 1, element = BiochemicalAlgorithms.Elements.O, name = "", atom_type = "", r = Float32[0.0, 0.0, 0.0], v = Float32[0.0, 0.0, 0.0], F = Float32[0.0, 0.0, 0.0], formal_charge = 0, charge = 0.0f0, radius = 0.0f0)
    Atom{Float32}: (idx = 3, number = 2, element = BiochemicalAlgorithms.Elements.H, name = "", atom_type = "", r = Float32[1.0, 0.0, 0.0], v = Float32[0.0, 0.0, 0.0], F = Float32[0.0, 0.0, 0.0], formal_charge = 0, charge = 0.0f0, radius = 0.0f0)
    Atom{Float32}: (idx = 4, number = 3, element = BiochemicalAlgorithms.Elements.H, name = "", atom_type = "", r = Float32[-0.25881904, 0.9659258, 0.0], v = Float32[0.0, 0.0, 0.0], F = Float32[0.0, 0.0, 0.0], formal_charge = 0, charge = 0.0f0, radius = 0.0f0)

    3-element Vector{Nothing}:
     nothing
     nothing
     nothing

# How can I determine the element type of an atom (i.e., C,N,…)?

The element of an atom is represented by \[`elements.jl`\]@ref

``` julia
s = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

for atom in atoms(s)
    if atom.element == Elements.C
        println("This is a carbon atom!")
    end
end
```

    This is a carbon atom!
    This is a carbon atom!
    This is a carbon atom!
    This is a carbon atom!
    This is a carbon atom!
    This is a carbon atom!

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

    N
    This is a backbone atom!
    C
    C
    This is a backbone atom!
    O
    This is a backbone atom!
    H
    H
    H
    H
    H
    H
    C
    H
    This is a backbone atom!
    N
    This is a backbone atom!
    C
    C
    This is a backbone atom!
    H
    H
    H
    H
    C
    H
    This is a backbone atom!
    O
    This is a backbone atom!
    O

# How can I create a peptide from its amino acid sequence?

# How can I determine the protein’s amino acids sequence in one or three-letter code?

# How can I get the one-letter code out of a pdb file?

s = load_pdb(ball_data_path(“../test/data/AlaAla.pdb”)) for chain in
chains(sys) for res in residues(sys) if res.
