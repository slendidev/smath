#!/usr/bin/env bash

set -euo pipefail

html_dir="${1:-build/docs/html}"
index_file="${html_dir%/}/index.html"

if [[ ! -f "${index_file}" ]]; then
  echo "error: ${index_file} not found" >&2
  exit 1
fi

perl -0pi -e 's{<div><div class="header">\s*<div class="headertitle"><div class="title">.*?</div></div>\s*</div><!--header-->\s*}{}s' "${index_file}"
