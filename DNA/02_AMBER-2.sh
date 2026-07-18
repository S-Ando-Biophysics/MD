#!/usr/bin/env bash
set -euo pipefail

cd 01_AMBER

cp ../leap.in .
tleap -f leap.in

cd ..