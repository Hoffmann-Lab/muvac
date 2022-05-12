#! /usr/bin/env bash
# (c) Konstantin Riege

src="$(readlink -e "$(dirname "$0")")"
[[ $# -eq 0 ]] && {
	exec "$src/bashbone/setup.sh" -h
} || {
	exec "$src/bashbone/setup.sh" -s "$src/lib/compile.sh" "$@"
}
