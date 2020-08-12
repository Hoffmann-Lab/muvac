#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && {
	[[ "${BASH_SOURCE[0]}" == "${0}" ]] && {
		echo ":ERROR: script needs to be sourced" >&2
		echo ":ERROR: do source $(readlink -e $0)" >&2
		exit 1
	} || {
		[[ ! $OSTYPE =~ linux ]] && echo "unsupported operating system" || {
			if [[ ${BASH_VERSINFO[0]} -gt 4 ]] || [[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]]; then
				insdir_pipeline=$(dirname $(readlink -e ${BASH_SOURCE[0]}))
				insdir_tools_pipeline=$(dirname $insdir_pipeline)
				unset OPTIND activate_conda_pipeline
				while getopts :i:c: arg; do
					case $arg in
						i) insdir_pipeline="$OPTARG";;
						c) activate_conda_pipeline="$OPTARG";;
						:) echo ":ERROR: argument missing" >&2 ; return 1;;
					esac
				done
				source "$insdir_pipeline/bashbone/activate.sh" -i "$insdir_tools_pipeline" -c ${activate_conda_pipeline:-false} || return 1
				BASHBONEVERSION=$version && \
				IFS=$'\n'
				for f in "$insdir_pipeline/lib/"*.sh; do
					source "$f"
				done && {
					unset IFS
					INSDIR="$insdir_pipeline"
					bashbone(){
						declare -f | grep -P '::.+\(\)' | grep -vF -e compile:: -e helper::_ -e progress:: -e commander::_ -e pipeline:: | sort -V | sed -r 's/\s+\(\)\s+$//'
					}
				} || {
					echo ":ERROR: pipeline activation failed" >&2
					echo ":ERROR: unexpected error in source code - please contact developer" >&2
					return 1
				}
			else
				echo ":ERORR: requieres bash version 4.4 or above" >&2
				return 1
			fi
		}
	}
} || {
	echo ":ERORR: loading library requieres bash" >&2
	return 1
}
