# Benchmark Report for */lustre/project/ki-cbr1fep/juliacon/BiochemicalAlgorithms.jl*

## Job Properties
* Time of benchmark: 19 Dec 2024 - 9:14
* Package commit: dirty
* Julia commit: bd47ec
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                                      | time            | GC time    | memory          | allocations |
|-------------------------------------------------------------------------|----------------:|-----------:|----------------:|------------:|
| `["FileFormats", "Json", "Reading"]`                                    | 452.036 μs (5%) |            | 671.23 KiB (1%) |        6529 |
| `["FileFormats", "PDB", "Reading"]`                                     |   7.719 ms (5%) |            |   8.47 MiB (1%) |      114992 |
| `["FileFormats", "SDfile", "Reading"]`                                  |  11.799 ms (5%) |            |  13.56 MiB (1%) |      100594 |
| `["ForceFields", "AmberFF", "Creation"]`                                | 111.534 ms (5%) |            |  90.26 MiB (1%) |      403125 |
| `["ForceFields", "AmberFF", "compute_energy!(::ForceField)"]`           |   4.358 ms (5%) |            | 234.70 KiB (1%) |        7493 |
| `["ForceFields", "AmberFF", "compute_energy!(::NonBonded)"]`            |   3.564 ms (5%) |            |                 |             |
| `["ForceFields", "AmberFF", "compute_energy!(::QuadraticAngleBend)"]`   |  72.679 μs (5%) |            |                 |             |
| `["ForceFields", "AmberFF", "compute_energy!(::QuadraticBondStretch)"]` |  15.650 μs (5%) |            |                 |             |
| `["ForceFields", "AmberFF", "compute_energy!(::Torsion)"]`              | 736.203 μs (5%) |            | 233.08 KiB (1%) |        7431 |
| `["ForceFields", "AmberFF", "compute_forces!(::ForceField)"]`           |  32.909 ms (5%) |            |  23.26 MiB (1%) |       10775 |
| `["ForceFields", "AmberFF", "setup!"]`                                  |  40.473 ms (5%) |            |  21.04 MiB (1%) |      351300 |
| `["ForceFields", "AmberFF", "update!"]`                                 |  62.905 ms (5%) |            |  66.66 MiB (1%) |        3694 |
| `["FragmentDB", "Building_Bonds"]`                                      |   6.270 ms (5%) | 262.993 μs |   5.04 MiB (1%) |      149239 |
| `["FragmentDB", "Creation"]`                                            |  56.843 ms (5%) |            |   4.65 MiB (1%) |      109103 |
| `["FragmentDB", "Normalize_Name"]`                                      |   6.211 ms (5%) | 258.251 μs |   4.98 MiB (1%) |      147368 |
| `["FragmentDB", "Reconstruct_Fragments"]`                               |   4.269 ms (5%) |            |   4.21 MiB (1%) |       70605 |
| `["Kernel", "Clone", "SystemCloning with bonds"]`                       |   3.430 ms (5%) |            |   2.42 MiB (1%) |       20791 |
| `["Kernel", "Clone", "SystemCloning wo bonds"]`                         |   2.610 ms (5%) |            |   1.19 MiB (1%) |       10689 |
| `["Kernel", "Creation", "AtomCreation"]`                                | 287.191 ns (5%) |            |  896 bytes (1%) |           8 |
| `["Kernel", "Creation", "FragmentCreation"]`                            | 262.538 ns (5%) |            |  896 bytes (1%) |           8 |
| `["Kernel", "Creation", "MoleculeCreation"]`                            | 195.911 ns (5%) |            |  896 bytes (1%) |           8 |
| `["Kernel", "Creation", "ResidueCreation"]`                             | 351.705 ns (5%) |            |  896 bytes (1%) |           8 |
| `["Kernel", "Creation", "SystemCreation"]`                              |   1.731 μs (5%) |            |   7.48 KiB (1%) |          83 |
| `["Kernel", "Iteration", "AtomIteration"]`                              |  42.340 μs (5%) |            |  34.89 KiB (1%) |         840 |
| `["Kernel", "Iteration", "BondIteration"]`                              | 112.279 μs (5%) |            |  95.88 KiB (1%) |        2167 |
| `["Kernel", "Iteration", "ChainIteration"]`                             |  20.611 ns (5%) |            |   64 bytes (1%) |           1 |
| `["Kernel", "Iteration", "ResidueIteration"]`                           |   2.319 μs (5%) |            |   1.92 KiB (1%) |           4 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["FileFormats", "Json"]`
- `["FileFormats", "PDB"]`
- `["FileFormats", "SDfile"]`
- `["ForceFields", "AmberFF"]`
- `["FragmentDB"]`
- `["Kernel", "Clone"]`
- `["Kernel", "Creation"]`
- `["Kernel", "Iteration"]`

## Julia versioninfo
```
Julia Version 1.10.2
Commit bd47eca2c8a (2024-03-01 10:14 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  uname: Linux 4.18.0-553.27.1.el8_10.x86_64 #1 SMP Tue Nov 5 04:50:16 EST 2024 x86_64 x86_64
  CPU: AMD EPYC 7713 64-Core Processor: 
                  speed         user         nice          sys         idle          irq
       #1-128  2000 MHz  2434117206 s     280534 s   12294388 s  522156874 s   14665768 s
  Memory: 251.66006088256836 GB (235109.08984375 MB free)
  Uptime: 2.33238201e6 sec
  Load Avg:  1.08  1.22  2.5
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, znver3)
Threads: 1 default, 0 interactive, 1 GC (on 128 virtual cores)
```