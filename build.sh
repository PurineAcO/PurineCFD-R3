#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

LOG="$(mktemp)"
if ! { cmake --preset clang-debug && cmake --build --preset clang-debug; } >"$LOG" 2>&1; then
    cat "$LOG" >&2
    rm -f "$LOG"
    exit 1
fi
rm -f "$LOG"

"./out/build/clang-debug/PurineCFD-R3.exe"
