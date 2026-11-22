# Maintainer checks

Run these commands from the package root. The R checks use Air 0.13.0 and lintr 3.4.0.
Air handles layout. Lintr checks correctness, allowing qualified Matrix calls and duplicate metric
names in lists. R CMD check supplies package-level usage diagnostics.

```sh
air format . --check
Rscript --vanilla tools/lint.R
python3 ~/.agents/skills/r-docs-pkgdown/scripts/format-rmd.py --check vignettes
```

Use `tools/native-quality.sh` on Linux with the package's LinkingTo dependencies installed,
GCC, Clang with libc++, clang-format 21.1.8 and clang-tidy. It compiles every installed header
independently and checks handwritten translation units, then exercises worker exceptions with
each C++ runtime and with both loaded. Set `CLANGXX`, `CLANG_FORMAT` or `CLANG_TIDY` to select
versioned executables. CI uses LLVM 21.

## Valgrind

Build in a scratch directory and install the resulting tarball in a fresh library. This makes
the smoke test print and verify the path of the package it loads. The check covers an Annoy fit,
model persistence, projection and unloading. Run the full package check for examples, tests
and vignette code.

```sh
repo=$(pwd)
work=$(mktemp -d)
mkdir "$work/lib"
cd "$work"
R CMD build "$repo"
R CMD INSTALL --library="$work/lib" uwot_*.tar.gz
R_LIBS="$work/lib" R -d 'valgrind --tool=memcheck --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=99' \
  --vanilla --slave -f "$repo/tools/valgrind-smoke.R" --args "$work/lib"
R CMD check --no-manual --use-valgrind uwot_*.tar.gz
```

Read Valgrind's error summaries as well as the R check result. The GitHub Actions R-hub workflow
runs weekly Clang and GCC checks and accepts the broader compiler, sanitizer, Valgrind, LTO and
rchk matrix on manual dispatch.
