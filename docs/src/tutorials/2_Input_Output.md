# All on Input and Output


BiochemicalAlgorithms.jl supports the reading and writing of several
common structural data formats.

## Protein Data Bank Format (PDB)

The most common format ist the PDB. Have a look at [Learning about PDB
data](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction)
if you want to refresh your knowledge about this format. For a deeper
understanding you can also read [PDB format
specifications](https://mmcif.wwpdb.org/pdbx-mmcif-home-page.html).

``` julia
using BiochemicalAlgorithms
```

``` julia
# Read PDB file from the BiochemicalAlgorithms.jl repository
pdb_sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))

for atom in atoms(pdb_sys)
    if atom.element == Elements.S
        print(atom.number)
    end
end

fdb = FragmentDB()

normalize_names!(pdb_sys,fdb)
reconstruct_fragments!(pdb_sys, fdb)
build_bonds!(pdb_sys, fdb)
write_pdb(ball_data_path("../test/data/Ala_out.pdb"), pdb_sys)
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 22 bonds

## MMCIF

We can also read and write mmcif files:

In addition to PDB, the [pubchem data
base](https://pubchem.ncbi.nlm.nih.gov/) plays an important role as a
source of structural data. Pubchem allows to retrieve data in JSON which
is read by BiochemicalAlgorithms.jl as shown below:

``` julia
pubchem_sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))

fdb = FragmentDB()
normalize_names!(pdb_sys,fdb)
reconstruct_fragments!(pdb_sys, fdb)
build_bonds!(pdb_sys, fdb)

for atom in atoms(pdb_sys)
    if atom.element == Elements.S
        print(atom.number)
    end
end
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 0 bonds

## SD file

``` julia
sd_sys = load_sdfile(ball_data_path("../test/data/rings_test.sdf"))

fdb = FragmentDB()
normalize_names!(sd_sys,fdb)
reconstruct_fragments!(sd_sys, fdb)
build_bonds!(sd_sys, fdb)

for atom in atoms(sd_sys)
    if atom.element == Elements.S
        print(atom.number)
    end
end

write_sdfile(ball_data_path("../test/data/rings_test_out.sdf"), sd_sys)
```

    [ Info: reconstruct_fragments!(): added 0 atoms.
    [ Info: build_bonds!(): built 0 bonds
    ┌ Warning: write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane...
    └ @ BiochemicalAlgorithms C:\local\BiochemicalAlgorithms.jl\src\fileformats\sdfile.jl:23
    [ Info: 9 records exported.
