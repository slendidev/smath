#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

doxygen "${repo_root}/docs/Doxyfile"
"${repo_root}/docs/postprocess-doxygen.sh" "${repo_root}/build/docs/html"

echo "Docs generated at ${repo_root}/build/docs/html"
