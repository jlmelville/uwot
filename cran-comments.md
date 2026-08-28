## CRAN-requested update

This release fixes issues reported for uwot 0.2.4 in the CRAN checks and the
additional Clang 23 checks.

## Test environments

* local Ubuntu 26.04 LTS, R 4.5.2, GCC 15.2.0
* local Ubuntu 26.04 LTS, R 4.5.2, Clang 21.1.8 with libstdc++ and libc++
* Ubuntu 24.04 (on GitHub Actions), R 4.5.3, R 4.6.1, and R-devel
  (4.7.0), GCC 13.3.0
* local Windows 11 build, R 4.6.1
* Windows Server 2025 with Visual Studio 2026 (on GitHub Actions), R 4.5.3
  and R 4.6.1
* Windows Server 2022 (on Win-Builder), R-devel (2026-08-27 r90452), GCC 14.3.0
* local macOS Tahoe R 4.6.1
* macOS Tahoe 26, arm64 (on GitHub Actions), R 4.6.1, Apple Clang 21.0.0
* Ubuntu 22.04 (on R-hub), R-devel 4.7.0, Clang 15.0.7
* Ubuntu 22.04 (on R-hub), R-devel 4.7.0, Clang 22.1.8 with libc++,
  ASAN and UBSAN
* Fedora 42 (on R-hub), R-devel 4.7.0, GCC 15.2.1, ASAN and Valgrind
* Fedora 44 (on R-hub), R-devel 4.7.0, GCC 16.2.1
* Ubuntu 24.04 (on R-hub), R 4.6.1, GCC 13.3.0, LTO

## R CMD check results

0 errors | 0 warnings | 0 notes

## CRAN checks

The CRAN checks for uwot 0.2.4 report one NOTE on the five r-devel flavors:

```
Version: 0.2.4
Check: R code for possible problems
Result: NOTE
  Found calls to structure() using deprecated special names:
    uwot/tests/testthat/test_output.R (.Dim: 1)
  '.Dim' should be changed to 'dim'.
```

Flavors: r-devel-linux-x86_64-debian-clang,
r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang,
r-devel-linux-x86_64-fedora-gcc, r-devel-windows-x86_64.

The test fixture now uses `dim` rather than `.Dim`.

The additional R-devel check on Fedora 44 with Clang 23.1.0 and libc++ reports
an installation ERROR:

```
Version: 0.2.4
Check: whether package can be installed
Result: ERROR
  ../inst/include/uwot/update.h:210:10: error: no member named 'fill' in namespace 'std'
  ../inst/include/RcppPerpendicular.h:97:52: error: no member named 'ref' in namespace 'std'
```

The installed C++ headers now directly include `<algorithm>` and `<functional>`,
which declare `std::fill` and `std::ref`, respectively. This avoids relying on
transitive libc++ includes removed in LLVM 23.

## revdepcheck results

We checked 98 reverse dependencies (58 from CRAN + 40 from Bioconductor),
comparing R CMD check results across CRAN and development versions of uwot.

* We saw 0 new problems.
* We failed to check 0 packages.
