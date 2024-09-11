using BenchmarkTools
using BiochemicalAlgorithms

const SUITE = BenchmarkGroup()

include("AmberFF_bench.jl")
include("KernelIteration_bench.jl")
include("KernelCreation_bench.jl")
include("FileFormats_bench.jl")
include("FragmentDB_bench.jl")
