# Benchmark Report for *..*

## Job Properties
* Time of benchmarks:
    - Target: 6 Feb 2024 - 14:25
    - Baseline: 6 Feb 2024 - 14:27
* Package commits:
    - Target: e9353d
    - Baseline: c5d009
* Julia commits:
    - Target: 312098
    - Baseline: 312098
* Julia command flags:
    - Target: None
    - Baseline: None
* Environment variables:
    - Target: None
    - Baseline: None

## Results
A ratio greater than `1.0` denotes a possible regression (marked with :x:), while a ratio less
than `1.0` denotes a possible improvement (marked with :white_check_mark:). Only significant results - results
that indicate possible regressions or improvements - are shown below (thus, an empty table means that all
benchmark results remained invariant between builds).

| ID                                             | time ratio                   | memory ratio                 |
|------------------------------------------------|------------------------------|------------------------------|
| `["ForceFields", "AmberFF", "Creation"]`       | 0.45 (5%) :white_check_mark: |                1.79 (1%) :x: |
| `["ForceFields", "AmberFF", "compute_energy"]` |                3.45 (5%) :x: |                   1.00 (1%)  |
| `["ForceFields", "AmberFF", "compute_forces"]` | 0.03 (5%) :white_check_mark: | 0.18 (1%) :white_check_mark: |
| `["ForceFields", "AmberFF", "setup"]`          | 0.14 (5%) :white_check_mark: | 0.40 (1%) :white_check_mark: |
| `["ForceFields", "AmberFF", "update!"]`        | 0.74 (5%) :white_check_mark: |                2.64 (1%) :x: |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["ForceFields", "AmberFF"]`

## Julia versioninfo

### Target
```
Julia Version 1.10.0
Commit 3120989f39b (2023-12-25 18:01 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin22.4.0)
  uname: Darwin 23.3.0 Darwin Kernel Version 23.3.0: Wed Dec 20 21:33:31 PST 2023; root:xnu-10002.81.5~7/RELEASE_ARM64_T8112 arm64 arm
  CPU: Apple M2: 
              speed         user         nice          sys         idle          irq
       #1  2400 MHz     694058 s          0 s     270616 s    1131016 s          0 s
       #2  2400 MHz     671612 s          0 s     268509 s    1174228 s          0 s
       #3  2400 MHz     625650 s          0 s     263218 s    1235854 s          0 s
       #4  2400 MHz     592397 s          0 s     250178 s    1290369 s          0 s
       #5  2400 MHz     173187 s          0 s     117132 s    1872489 s          0 s
       #6  2400 MHz     135094 s          0 s      72103 s    1967074 s          0 s
       #7  2400 MHz      94988 s          0 s      46130 s    2037807 s          0 s
       #8  2400 MHz      81882 s          0 s      33993 s    2065719 s          0 s
  Memory: 24.0 GB (111.5625 MB free)
  Uptime: 269305.0 sec
  Load Avg:  14.8974609375  13.11279296875  13.39990234375
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, apple-m1)
  Threads: 1 on 4 virtual cores
```

### Baseline
```
Julia Version 1.10.0
Commit 3120989f39b (2023-12-25 18:01 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin22.4.0)
  uname: Darwin 23.3.0 Darwin Kernel Version 23.3.0: Wed Dec 20 21:33:31 PST 2023; root:xnu-10002.81.5~7/RELEASE_ARM64_T8112 arm64 arm
  CPU: Apple M2: 
              speed         user         nice          sys         idle          irq
       #1  2400 MHz     694899 s          0 s     270986 s    1131054 s          0 s
       #2  2400 MHz     672469 s          0 s     268910 s    1174266 s          0 s
       #3  2400 MHz     626506 s          0 s     263642 s    1235890 s          0 s
       #4  2400 MHz     593253 s          0 s     250611 s    1290406 s          0 s
       #5  2400 MHz     174015 s          0 s     117570 s    1872561 s          0 s
       #6  2400 MHz     135904 s          0 s      72543 s    1967144 s          0 s
       #7  2400 MHz      95806 s          0 s      46571 s    2037884 s          0 s
       #8  2400 MHz      82693 s          0 s      34435 s    2065802 s          0 s
  Memory: 24.0 GB (86.03125 MB free)
  Uptime: 269446.0 sec
  Load Avg:  47.955078125  27.642578125  19.21875
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, apple-m1)
  Threads: 1 on 4 virtual cores
```