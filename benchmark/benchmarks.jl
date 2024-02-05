using BenchmarkTools
using BiochemicalAlgorithms

const SUITE = BenchmarkGroup()

include("AmberFF_bench.jl")