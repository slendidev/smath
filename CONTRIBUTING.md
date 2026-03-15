# Contributing to smath

Welcome to `smath` contributors!

This document contains a set of guidelines to contribute to the project. If you feel like there's a mistake or something can be improved to this document, feel free to propose changes in a pull request.

## Philosophy

- `smath` itself should be kept header only, in a single file.
- `smath` should be simple and easy to use.

## Development setup

You can work with either Nix (recommended, same tooling used in CI) or a local CMake toolchain.

### Option 1: Nix

Make sure you have flakes enabled and then execute the following to enter the shell:

```bash
nix develop
```

You can also run `direnv allow` if you have [direnv](https://direnv.net/) installed and configured.

### Option 2: Local toolchain

Requirements:

- C++23-capable compiler
- CMake (3.15+)
- Ninja (recommended)

### Configure and build

```bash
cmake -S . -B build -G Ninja -DSMATH_BUILD_TESTS=ON -DSMATH_BUILD_EXAMPLES=ON
cmake --build build
```

## Running tests

```bash
ctest --test-dir build --output-on-failure
```

## Code style

This project uses `clang-format` to maintain a consistent code style. Before opening a pull request, format changed C++ files:

```bash
git ls-files '*.hpp' '*.cpp' | xargs clang-format -i
```

## Pull requests

- Keep changes focused and scoped to one topic.
- Add or update tests when behavior changes.
- Ensure the project builds and tests pass locally.

### Commit style

Use the commit message format: `<category>: <brief>`.

Allowed categories: `feat`, `fix`, `test`, `docs`, `ci`.

Examples:

- `feat: add as_matrix() to quaternion`
- `fix: correct mat4 approx_equal ignoring a column`
- `test: add tests for vector swizzle edge cases`
- `docs: fix typo in README.md`
- `ci: update nix build command in pull request workflow`

If a commit fixes a tracked issue, include an issue-closing footer: `Closes: #<num>` (example: `Closes: #42`).

## Reporting issues

When opening an issue, include:

- what you expected to happen
- what happened instead
- a minimal reproducible example
- compiler and platform details
