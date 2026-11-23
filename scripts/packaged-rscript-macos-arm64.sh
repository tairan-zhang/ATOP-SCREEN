#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
software_root="$(cd "${script_dir}/.." && pwd)"
r_home="${software_root}/packaging/runtime/macos-arm64/runtime/R.framework/Versions/4.6/Resources"
r_exec="${r_home}/bin/exec/R"

if [[ ! -x "${r_exec}" ]]; then
  echo "Missing packaged R executable: ${r_exec}" >&2
  exit 1
fi
if [[ $# -lt 1 ]]; then
  echo "Usage: packaged-rscript-macos-arm64.sh script.R [arguments ...]" >&2
  exit 1
fi

script_path="$1"
shift

export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"
export R_HOME="${r_home}"
if [[ -n "${ATOP_EXTRA_R_LIBS:-}" ]]; then
  export R_LIBS="${ATOP_EXTRA_R_LIBS}:${r_home}/library"
  export R_LIBS_USER="${ATOP_EXTRA_R_LIBS}:${r_home}/library"
else
  export R_LIBS="${r_home}/library"
  export R_LIBS_USER="${r_home}/library"
fi
export R_SHARE_DIR="${r_home}/share"
export R_INCLUDE_DIR="${r_home}/include"
export R_DOC_DIR="${r_home}/doc"
export PATH="${r_home}/bin:${PATH}"
export DYLD_LIBRARY_PATH="${r_home}/lib:${DYLD_LIBRARY_PATH:-}"

exec "${r_exec}" --no-echo --no-restore --vanilla --file="${script_path}" --args "$@"
