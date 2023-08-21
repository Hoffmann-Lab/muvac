#! /usr/bin/env bash
# (c) Konstantin Riege

src="$(dirname "$(dirname "$(readlink -e "$0")")")"
[[ $# -eq 0 ]] && {
	exec "$src/bashbone/scripts/setup.sh" -h
} || {
	exec "$src/bashbone/scripts/setup.sh" -s "$src/lib/compile.sh" "$@"
}
