# All on Input and Output


BiochemicalAlgorithms.jl supports the reading and writing of several common structural data formats.

## Protein Data Bank Format (PDB)

The most common format ist the PDB. Have a look at [Learning about PDB data](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction) if you want to refresh your knowledge about this format. For a deeper understanding you can also read [PDB format specifications](https://mmcif.wwpdb.org/pdbx-mmcif-home-page.html).

Read a PDB file from the BiochemicalAlgorithms.jl repository:

``` julia
sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"))
```

    System{Float32}: AlaAla.pdb
      23 atoms
       0 bonds
       1 molecules
       1 chains
       0 secondary structures
       2 fragments

Write the same system back into a new PDB file:

``` julia
write_pdb("Ala_out.pdb", sys)
```

## PubChem

In addition to PDB, the [pubchem data base](https://pubchem.ncbi.nlm.nih.gov/) plays an important role as a source of structural data. Pubchem allows to retrieve data in JSON which is read by BiochemicalAlgorithms.jl as shown below:

``` julia
sys = load_pubchem_json(ball_data_path("../test/data/aspirin_pug.json"))
```

    System{Float32}:
      21 atoms
      21 bonds
       1 molecules
       0 chains
       0 secondary structures
       0 fragments

## SD file

``` julia
sys = load_sdfile(ball_data_path("../test/data/rings_test.sdf"))
```

    System{Float32}: rings_test.sdf
      148 atoms
      151 bonds
        9 molecules
        0 chains
        0 secondary structures
        0 fragments

Write the system into a new SD file:

``` julia
write_sdfile("rings_test_out.sdf", sys)
```

    ┌ Warning: write_sdfile: writer only supports 2D data; projecting atoms onto xy-plane...
    └ @ BiochemicalAlgorithms ~/git/ball.jl/src/fileformats/sdfile.jl:33
    [ Info: 9 records exported.
