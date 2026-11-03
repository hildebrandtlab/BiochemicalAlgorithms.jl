# Information for AI agents

## Module Organization
- [`src/core/`](src/core/) - `System` and its components (atoms, bonds, molecules, chains, fragments, secondary structures), plus their tables
- [`src/fileformats/`](src/fileformats/) - Readers and writers for PDB, mmCIF, HIN, SDF, PubChem JSON, and BALL INI files
- [`src/forcefields/`](src/forcefields/) - Shared force field components (`common/`), AMBER, and MMFF94 parameters
- [`src/mappings/`](src/mappings/) - Atom bijections and rigid mappings
- [`src/optimization/`](src/optimization/) - Structure optimization
- [`src/preprocessing/`](src/preprocessing/) - `FragmentDB`, name normalization, fragment reconstruction, bond building, and hydrogen bond and secondary structure prediction
- [`src/structureanalysis/`](src/structureanalysis/) - Bounding boxes and molecular surfaces
- [`src/substructures/`](src/substructures/) - Substructures, SMARTS, and SSSR
- [`data/`](data/) - Force field and fragment data shipped with the package (see `ball_data_path`)
- [`test/`](test/) - Test items, laid out like `src/`, and test structures in `test/data/`
- [`docs/`](docs/) - Documenter setup and manual pages
- [`tutorials/`](tutorials/) - Quarto sources of the tutorials
- [`benchmark/`](benchmark/) - PkgBenchmark suite

All source files are included from [`src/BiochemicalAlgorithms.jl`](src/BiochemicalAlgorithms.jl). The order
matters, so add new files at the point where their dependencies are already defined.

## Running Julia

You should have a recent binary copy of julia in your `$HOME/.juliaup/bin` directory.
You may use this julia executable for validation. The package supports Julia 1.10 and later, so do
not use language or library features that are newer than 1.10, even if your local julia has them.

The first precompilation of the package takes about a minute. To avoid repeating it while iterating,
load `Revise` before the package in a long-running session:

```julia
using Revise, BiochemicalAlgorithms
```

## For all changes

1. Sign off every commit you create with an `Assisted-by: <tool> (<model>)` trailer naming both the
   tool and the model behind it, for example `Assisted-by: Claude Code (Opus 5)`. Name the model,
   not just the harness - it is what tells a later reader what actually produced the work. Do not
   use `Co-authored-by:` for tools, and do not list a tool as an author; this trailer replaces any
   co-author trailer your harness adds by default.
2. Every commit also needs a `Signed-off-by: Name <E-mail>` line, which CI checks. This line certifies
   the [Developer Certificate of Origin](CONTRIBUTING.md) on behalf of a human, so only add it (for
   example with `git commit -s`) when the human user asks you to, and never in the name of a tool.
3. Agents can only open PRs or post comments once the human user gives them explicit permission.
4. Do not commit `Manifest.toml` files; they are ignored on purpose.

## Testing

Tests are written as `@testitem`s and discovered automatically from any `.jl` file under `test/`, so
they do not need to be registered anywhere. To run the whole test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test dependencies are not part of a separate `test/Project.toml`. To run a subset of the tests, use
a temporary environment from the repository root and filter by `filename`, `name`, or `tags`:

```bash
julia -e '
using Pkg
Pkg.activate(temp=true)
Pkg.develop(path=".")
Pkg.add(["TestItemRunner", "Aqua", "Test"])
using TestItemRunner
TestItemRunner.run_tests("."; filter=ti->occursin("test_bond_order", ti.filename), verbose=true)
'
```

When writing tests:

- Each `@testitem` runs in its own module with `BiochemicalAlgorithms` and `Test` already loaded. Any
  other package has to be loaded inside of the test item.
- Cover both `Float32` and `Float64` (`for T in [Float32, Float64]`) when the code is generic over the
  floating-point type.
- Load test structures with `ball_data_path("../test/data/<file>")`.
- Tag test items with `:skip_ci` if they should not run on GitHub Actions.
- Do not change expected energies or forces just to make a test pass. Many of these values are reference
  values, some taken from the [BALL](https://github.com/BALL-Project/ball) test suite. If a change alters
  them, explain why to the human user.

[`test/aqua.jl`](test/aqua.jl) runs the Aqua checks as part of the test suite. Among other things, every new
dependency needs a `[compat]` entry in `Project.toml`.

## Building the documentation

Build the manual locally with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

- Docstrings have to be listed in an `@docs` block to end up in the manual: exported names in
  `docs/src/public/`, and internal names in `docs/src/private/`.
- When changing a tutorial, edit the `.qmd` source in `tutorials/`. The Markdown files in
  `docs/src/tutorials/` are generated from it, but are also tracked, so regenerate them by running
  `quarto render` in `tutorials/` and commit both. Do not edit only the generated file.

## Running benchmarks

If your change affects performance-sensitive code, such as force field evaluation or optimization, compare
the benchmarks before and after the change. Run them from the `benchmark/` directory:

```bash
cd benchmark && julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate(); include("run_benchmarks.jl")'
```

Results are written to `benchmark_results.md`, which is ignored by git.

## Writing code

After writing code, look up the docstring for each function you used. If there
are recommendations or additional considerations that apply to these functions,
make sure to take them into account.

- Keep code generic over the floating-point type (`T <: Real`) rather than hardcoding `Float32` or `Float64`.
  Note that `System()` defaults to `Float32`.
- List exported names in an `export` block at the top of the file that defines them.
- Prefix internal functions with an underscore (e.g. `_next_idx!`).
- Write docstrings for new API, and use `$(TYPEDSIGNATURES)` and `$(TYPEDEF)` from DocStringExtensions
  where the surrounding code does.
- Do not explicitly add a module prefix if the code you're adding is in the same module. E.g. do not use
  `BiochemicalAlgorithms.` for code in BiochemicalAlgorithms unless required.

## Commit messages and pull requests

When writing commit messages, follow the format "component: Brief summary" for
the title, where the component is, for example, `PDB`, `mmCIF`, `System`, `FragmentDB`, `Docs`, or `CI`.
In the body of the commit message, provide a brief prose summary
of the purpose of the changes made. Do not specifically mention added tests, comments,
documentation, etc., unless this is the main purpose of the change. Do not mention
the test plan, unless it differs from what you were instructed to do in AGENTS.md.
If your change fixes one or more issues, use the syntax "Fixes #" at the end of the commit message, but do not include it in the title.

When referencing external GitHub PRs or issues, use proper GitHub interlinking format (e.g., `owner/repo#123` for PRs/issues).
When fixing CI failures, include the link to the specific CI failure in the commit message.
Always quote macro names in backticks in commit messages and PR titles/bodies (e.g. `` `@testitem` ``, not @testitem), so GitHub does not notify the unrelated user with that handle.

When preparing a pull request for the human author to open, draft the body as follows, expecting
them to reword it:
1. If the pull request consists of one commit only, use the body of the commit for the body of the pull request.
2. If there are multiple commits in the pull request, follow the same guidelines for the pull request as for the commit body.
3. Target the `develop` branch, and make sure that the pull request is based on a recent commit of it - if not rebase your changes first.
4. If a separate tool reviews the commit, this is useful to note. Encourage the human to state how carefully they read and understood the content also, for example, by drafting the text to say the human author has *not* read any of it, and expecting they will update that once they read that disclaimer.
