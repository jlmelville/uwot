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
* local Ubuntu 26.04 LTS, R 4.5.2, Clang 21.1.8 with libstdc++

External CRAN-like compiler, sanitizer, and platform checks are still pending.
Their actual environments will be added after the corresponding runs complete.

## R CMD check results

0 errors | 0 warnings | 0 notes

This result is from the local GCC environment above with `--as-cran` and
`--no-manual`. A separate source installation with Clang, `-Wall`, `-Wextra`,
and `-pedantic` also completed successfully.

## revdepcheck results

Pending a fresh reverse-dependency check against the release candidate.
