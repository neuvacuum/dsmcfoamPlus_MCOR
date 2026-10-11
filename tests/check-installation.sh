#!/usr/bin/env bash
set -eo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$root/etc/bashrc"
for app in dsmcFoam+ dsmcInitialise+; do
    binary=$FOAM_USER_APPBIN/$app
    test -x "$binary"
    links=$(ldd "$binary")
    printf '%s\n%s\n' "$binary" "$links"
    if grep -q 'not found' <<< "$links"; then
        echo "Unresolved runtime library." >&2
        exit 1
    fi
    for lib in liblagrangian+.so libgeneralMolecule.so libdsmcFoam+.so; do
        grep -F "$FOAM_USER_LIBBIN/$lib " <<< "$links" >/dev/null || {
            echo "Wrong or missing extension library: $lib" >&2
            exit 1
        }
    done
    "$binary" -help >/dev/null 2>&1
done
echo "MCOR_INSTALLATION_PASS"
