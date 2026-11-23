#!/usr/bin/env bash
# Run from the source root with conda and conda-pack 0.8.1 available.
set -euo pipefail
[[ "$(uname -s)-$(uname -m)" == Darwin-arm64 ]]
[[ "$(conda-pack --version)" == *0.8.1 ]]
target=packaging/runtime/macos-arm64/mageck
if [[ -e "$target" ]]; then
  echo "Runtime already exists at $target" >&2
  exit 1
fi
staging=$(mktemp -d)
trap 'rm -rf "$staging"' EXIT
conda create --yes --prefix "$staging/env" --file packaging/locks/mageck-osx-arm64.txt
conda-pack --prefix "$staging/env" --output "$staging/runtime.tar" --format tar
mkdir "$staging/bundle"
tar -xf "$staging/runtime.tar" -C "$staging/bundle"
cp packaging/sources/MAGeCK-LICENSE "$staging/bundle/"
[[ "$("$staging/bundle/bin/python" -IB "$staging/bundle/bin/mageck" --version)" == 0.5.9.5 ]]
mkdir -p "$(dirname "$target")"
mv "$staging/bundle" "$target"
