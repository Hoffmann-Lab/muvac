#! /usr/bin/env bash
# (c) Konstantin Riege

callvariants::vcfzip() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-i <vcf>      | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads vcf
	while getopts 'S:s:t:i:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			i) ((mandatory++)); vcf="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage && return 1

	commander::print "compressing vcf"

	declare -a cmd1 cmd2
	readlink -e "$vcf" | file -f - | grep -Eo 'gzip' &&	commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
		bgzip -f -@ $threads < "$vcf" > "$vcf.gz"
	CMD
	commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
		tabix -f -p vcf "$vcf.gz"
	CMD

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t $threads -a cmd1
			commander::runcmd -v -b -t $threads -a cmd1
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
