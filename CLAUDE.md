# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A pure C99 numerical computing library for dynamics problems, focusing on Kalman filtering and linear algebra operations. Currently implementing an Unscented Kalman Filter based on the Wan & Merwe paper.

## Build Commands

```bash
# Build
mkdir build && cd build
cmake ..
make

# Run tests
./tests
```

Build uses GCC with strict flags: `-std=gnu99 -Wall -Wextra -Wfatal-errors -Werror`

## Architecture

### Module Structure

The library is organized into four main modules, each with headers in `inc/` and implementations in `src/`:

- **data_structures/** - Core containers (`v_t` for vectors, `m_t` for matrices)
- **linear_algebra/** - Decompositions (QR, Cholesky, Gram-Schmidt) and property checks
- **filtering/** - Unscented Kalman Filter implementation
- **integrators/** - ODE integration interface (skeleton only)

### Dependencies Between Modules

```
vector, matrix (foundational, no internal deps)
    ↓
linear_algebra_properties
    ↓
linear_algebra_decompositions
    ↓
filtering_kalman
```

### Error Handling Pattern

Functions return `error_t` codes defined in `inc/errors.h`:
- `E_OK` (0) - Success
- `E_NULLP` (1) - Null pointer
- `E_VAL` (2) - Invalid value/dimension mismatch
- `E_ERR` (3) - Generic error
- `E_INIT` (4) - Initialization error

For floating-point returns, `V_NAN` signals errors.

### Naming Conventions

- Module prefixes: `v_` (vector), `m_` (matrix), `la_` (linear algebra), `kalman_`
- Types use `_t` suffix: `v_t`, `m_t`, `error_t`

### Memory Management

- Manual allocation via `_new()` functions, deallocation via `_free()`/`_del()`
- Callers create output containers and are responsible for freeing them
- Functions return NULL on allocation failure

## Testing

Uses the Clar test framework (git submodule in `tests/clar`). Test files mirror source structure.

Key test macros: `cl_assert()`, `cl_assert_equal_i_()`

## Key Files

- `src/data_structures/matrix.c` - Core matrix operations (largest module)
- `src/filtering/kalman.c` - Unscented Kalman Filter implementation
- `src/linear_algebra/decompositions.c` - QR, Cholesky, Gram-Schmidt algorithms
- `inc/filtering/kalman.h` - Kalman filter context and callback types
