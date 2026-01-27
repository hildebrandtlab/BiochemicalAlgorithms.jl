# PDB Reader Test Findings

## Overview

Extensive testing of `load_pdb()` against 1000+ real PDB structures from RCSB reveals several failure modes and areas for improvement.

## Test Results Summary

| Metric | Value |
|--------|-------|
| Structures tested | 1239 unique PDB IDs |
| Categories | 30 (enzymes, membrane proteins, nucleic acids, etc.) |
| Success rate | **~72%** |
| Parse failures | ~28% |

## Test Details

- **Test file:** `test/fileformats/test_pdb_extensive.jl`
- **Structures tested:** 1239 unique PDB IDs across 30 categories
- **Success rate:** ~72% (stable across 600+ structures tested)
- **Cache location:** `~/.cache/BiochemicalAlgorithms_PDB_test/`

## Unrecognized PDB Record Types

The parser logs warnings for these unrecognized record types (in `src/fileformats/pdb/pdb_general.jl:60`):

| Record | Description | Priority |
|--------|-------------|----------|
| `NUMMDL` | Number of models in NMR ensemble | Low |
| `MDLTYP` | Model type annotation | Low |
| `DBREF1` | Extended database reference (part 1) | Medium |
| `DBREF2` | Extended database reference (part 2) | Medium |

**Recommendation:** Add parsing support for `DBREF1`/`DBREF2` which are used in modern PDB files for long database identifiers (UniProt, etc.).

## Parse Failure Categories

### 1. Missing Element Data (~20% of failures)
```
ErrorException("type Nothing has no field element")
```
**Cause:** Atom records where element cannot be determined from atom name.
**Affected structures:** Often older PDB files or those with non-standard atom names.
**Fix location:** `src/fileformats/pdb/pdb_atom.jl` - element inference logic.

### 2. Empty Fragment/Chain Errors (~60% of failures)
```
MethodError(Base.mapreduce_empty, (BiochemicalAlgorithms.PDBDetails.var"#27#29"(),
  Base.BottomRF{typeof(min)}(min), Fragment{Float32}), ...)
```
**Cause:** Calling `min()`/`max()` on empty collections when processing fragments.
**Affected structures:** PDB files with unusual chain/residue organization.
**Fix location:** Add empty collection checks before `minimum()`/`maximum()` calls.

### 3. Other Errors (~20% of failures)
- `BoundsError` - Index out of bounds during parsing
- `KeyError` - Missing lookup values
- Various edge cases

## Category-Specific Results (from 200-structure test)

| Category | Success Rate |
|----------|-------------|
| Nucleic acids | 96.7% |
| Landmark proteins | 85.0% |
| Edge cases | 82.8% |
| Diverse proteins | 76.0% |
| Recent structures (2020+) | 38.5% |

**Note:** Recent structures have lower success rates due to newer PDB format features (DBREF1/DBREF2, etc.).

## Recommended Fixes (Priority Order)

1. **High:** Add empty collection guards before `min()`/`max()` operations
2. **High:** Improve element inference for non-standard atom names
3. **Medium:** Add `DBREF1`/`DBREF2` record parsing
4. **Low:** Add `NUMMDL`/`MDLTYP` record parsing

## Running the Test

```bash
# Run extensive test (requires network, ~30-60 min first run)
julia --project -e 'using TestItemRunner; @run_package_tests filter=ti->(ti.name == "Extensive PDB Reader")'

# Subsequent runs are faster due to caching
```

## Test Configuration

The test uses a 60% success threshold to allow CI-like validation while acknowledging known issues. Adjust in `test_pdb_extensive.jl`:
```julia
@test success_rate >= 0.60  # Adjust threshold as fixes are implemented
```
