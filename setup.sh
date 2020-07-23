#! /usr/bin/env bash
# (c) Konstantin Riege

src=$(readlink -e $(dirname $0))
[[ $# -eq 0 ]] && {
	$src/bashbone/setup.sh -h
} || {
	$src/bashbone/setup.sh -s $src/lib/compile.sh "$@"
}
exit $?
