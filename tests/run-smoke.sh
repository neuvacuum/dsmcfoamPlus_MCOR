#!/usr/bin/env bash
set -eo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$root/etc/bashrc"
"$root/tests/check-installation.sh"
if [[ $# -gt 1 ]]; then
    echo "Usage: tests/run-smoke.sh [new-output-directory]" >&2
    exit 2
fi
if [[ $# == 1 ]]; then
    mkdir -- "$1"
    case_dir=$(cd -- "$1" && pwd -P)
else
    mkdir -p "$root/runs"
    case_dir=$(mktemp -d "$root/runs/smoke.XXXXXXXX")
fi
trap 'status=$?; printf "MCOR_SMOKE_FAILED: inspect %s/log.*\n" "$case_dir" >&2; exit "$status"' ERR
cp -R "$root/examples/closedBox/." "$case_dir/"
blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1
"$FOAM_USER_APPBIN/dsmcInitialise+" -case "$case_dir" > "$case_dir/log.dsmcInitialise" 2>&1
"$FOAM_USER_APPBIN/dsmcFoam+" -case "$case_dir" > "$case_dir/log.dsmcFoam" 2>&1
grep -q 'End main' "$case_dir/log.dsmcFoam"
grep -q 'Moment-conservative occupancy (Cartesian)' "$case_dir/log.dsmcFoam"
if grep -Ei 'FOAM FATAL|segmentation fault|nan' "$case_dir/log.dsmcFoam"; then
    exit 1
fi
awk '
    /Total Iterations =/ { iterations = $4 }
    /^[[:space:]]*Collisions[[:space:]]*=/ { collisions += $3 }
    /^[[:space:]]*Collision candidates[[:space:]]*=/ { candidates += $4 }
    /Moment-conservative occupancy/ {
        stages++
        for (i = 1; i <= NF; i++) {
            split($i, field, "=")
            if (field[1] == "deletedParcels") deleted += field[2]
            if (field[1] == "maxScaledMomentResidual" && field[2] + 0 > residual)
                residual = field[2] + 0
        }
    }
    END {
        printf "iterations=%d stages=%d candidates=%d collisions=%d deleted=%d maxScaledMomentResidual=%.9g\n", iterations, stages, candidates, collisions, deleted, residual
        if (iterations != 10 || stages != 20 || candidates < 1 || collisions < 1 || deleted < 1 || residual > 1e-9)
            exit 1
    }
' "$case_dir/log.dsmcFoam" | tee "$case_dir/smoke-summary.txt"
printf 'MCOR_SMOKE_PASS\nCase: %s\n' "$case_dir"
