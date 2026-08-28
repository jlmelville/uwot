## CRAN-requested update

This update addresses the problems reported on the CRAN check results page and
in the additional Clang 23 checks:

* A test fixture now uses `dim` rather than the deprecated special name `.Dim`
  in a call to `structure()`.
* Installed C++ headers now include the standard-library headers that declare
  the facilities they use. This avoids relying on transitive libc++ includes
  removed in LLVM 23.

## Test environments

* local Ubuntu 26.04 LTS, R 4.5.2, GCC 15.2.0
* local Ubuntu 26.04 LTS, R 4.5.2, Clang 21.1.8 with libstdc++ and libc++
* Ubuntu 24.04 (on GitHub Actions), R 4.5.3, R 4.6.1, and R-devel
  (4.7.0), GCC 13.3.0
* local Windows 11 build, R 4.6.1
* Windows Server 2025 with Visual Studio 2026 (on GitHub Actions), R 4.5.3
  and R 4.6.1
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

## revdepcheck results

We checked 98 reverse dependencies (58 from CRAN + 40 from Bioconductor),
comparing R CMD check results across CRAN and development versions of uwot.

* We saw 0 new problems.
* We failed to check 0 packages.
