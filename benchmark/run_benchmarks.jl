using PkgBenchmark

result = benchmarkpkg("..")

export_markdown("benchmark_results.md", result)