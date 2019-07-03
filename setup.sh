#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
shopt -s extglob

die() {
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

export SRC=$(readlink -e $(dirname $0))
[[ $# -eq 0 ]] && {
	$SRC/bashbone/setup.sh -h
} || {
	$SRC/bashbone/setup.sh -s $SRC/lib/compile.sh "$@"
}
exit 0
