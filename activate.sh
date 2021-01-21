#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ $OSTYPE != "linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

insdir="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
toolsdir="$(dirname "$insdir")"
unset OPTIND activate
while getopts ':i:c:x:' arg; do
	case $arg in
		i) toolsdir="$OPTARG";;
		c) activate="$OPTARG";;
		x) exitfun="$OPTARG";;
		:) echo "argument missing" >&2; return 1;;
	esac
done

source "$insdir/bashbone/activate.sh" -i "$toolsdir" -c ${activate:-false} -x "$exitfun" || exit 1

_IFS=$IFS
IFS=$'\n'
for f in "$insdir/lib/"*.sh; do
	BASHBONE_ERROR="file not found $f"
	[[ -s "$f" ]]
	source "$f"
done
IFS=$_IFS

unset BASHBONE_ERROR
return 0
