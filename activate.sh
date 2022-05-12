#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

insdir="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
toolsdir="$(dirname "$insdir")"
unset OPTIND activate exitfun help
while getopts ':i:c:x:a:h' arg; do
	case $arg in
		i) toolsdir="$OPTARG";;
		c) activate="$OPTARG";;
		x) exitfun="$OPTARG";;
		a) shift $((OPTIND-2)); break;;
		h) help=true;;
	esac
done

if [[ $BASHBONE_DIR ]]; then
	source "$BASHBONE_DIR/activate.sh" ${help:+-h} -i "${BASHBONE_TOOLSDIR:-$toolsdir}" -c ${activate:-false} -x "$exitfun" -a "$@" || return 1
else
	source "$insdir/bashbone/activate.sh" ${help:+-h} -i "$toolsdir" -c ${activate:-false} -x "$exitfun" -a "$@" || return 1
fi

[[ "$PATH" =~ ^$insdir: ]] || PATH="$insdir:$PATH"

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
