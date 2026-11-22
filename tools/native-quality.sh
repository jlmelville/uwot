#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
clangxx=${CLANGXX:-clang++}
formatter=${CLANG_FORMAT:-clang-format}
tidy=${CLANG_TIDY:-clang-tidy}
if [[ ! $("$formatter" --version) =~ version\ 21\.1\.8([^0-9.]|$) ]]; then
  echo 'Install clang-format 21.1.8 or point CLANG_FORMAT to it.' >&2
  exit 1
fi

# R and LinkingTo headers are dependencies. Apply our warning policy to uwot.
includes=$(Rscript --vanilla -e '
  deps <- strsplit(read.dcf("DESCRIPTION")[1, "LinkingTo"], ",")[[1]]
  deps <- trimws(sub(" *\\(.*", "", deps))
  paths <- vapply(deps, function(p) system.file("include", package = p,
                                                mustWork = TRUE), "")
  cat(c(R.home("include"), paths), sep = "\n")
')
flags=(-std=c++17 -pthread -fPIC -DRCPP_NO_RTTI -Isrc -Iinst/include)
while IFS= read -r path; do
  flags+=(-isystem "$path")
done <<< "$includes"
warnings=(-Wall -Wextra -Wpedantic -Wformat=2 -Wnull-dereference -Werror)

shopt -s globstar nullglob
headers=(inst/include/**/*.h)
sources=()
for source in src/*.cpp; do
  if [[ $source != src/RcppExports.cpp ]]; then
    sources+=("$source")
  fi
done
"$formatter" --dry-run --Werror "${sources[@]}" src/*.h "${headers[@]}" tools/*.cpp

scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
for compiler in g++ "$clangxx"; do
  library=()
  if [[ $compiler == "$clangxx" ]]; then
    library=(-stdlib=libc++)
  fi
  "$compiler" --version | head -n 1
  for header in "${headers[@]}"; do
    printf '#include <%s>\n' "${header#inst/include/}" |
      "$compiler" "${library[@]}" "${flags[@]}" "${warnings[@]}" -x c++ -fsyntax-only -
  done
  "$compiler" "${library[@]}" "${flags[@]}" "${warnings[@]}" -fsyntax-only "${sources[@]}"
  "$compiler" "${library[@]}" "${flags[@]}" "${warnings[@]}" \
    tools/test-worker-exceptions.cpp -o "$scratch/workers"
  "$scratch/workers"
done

# R can bring libstdc++ into a process that uses libc++ for the package.
libstdcpp=$(g++ -print-file-name=libstdc++.so.6)
test -f "$libstdcpp"
LD_PRELOAD="$libstdcpp" "$scratch/workers"

"$tidy" --warnings-as-errors='*' \
  --checks='-*,clang-analyzer-core.*,bugprone-unused-raii,bugprone-use-after-move,performance-unnecessary-copy-initialization' \
  --header-filter='(^|/)(src|inst/include)/' "${sources[@]}" -- "${flags[@]}"
