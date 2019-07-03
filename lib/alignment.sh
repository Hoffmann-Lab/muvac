#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::addreadgroup() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (requieres -2)
			-2 <tumoridx>  | array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir outdir i
	declare -n _mapper_addreadgroup _bamslices_addreadgroup _nidx_addreadgroup _tidx_addreadgroup
	declare -A nidx tidx
	while getopts 'S:s:t:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			r) ((mandatory++)); _mapper_addreadgroup=$OPTARG;;
			1) _nidx_addreadgroup=$OPTARG;;
			2) _tidx_addreadgroup=$OPTARG;;
			c) ((mandatory++)); _bamslices_addreadgroup=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1
	if [[ ! $_nidx_addreadgroup ]]; then
		declare -n _bams_addreadgroup="${_mapper_addreadgroup[0]}"
		for i in "${!_bams_addreadgroup[@]}"; do
			nidx[$i]=1
		done
	else
		[[ ! $_tidx_addreadgroup ]] && _usage && return 1
		for i in "${!_nidx_addreadgroup[@]}"; do
			nidx[$i]=1
		done
		for i in "${!_tidx_addreadgroup[@]}"; do
			tidx[$i]=1
		done
	fi

	commander::print "replacing read group tags"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o rgprefix slice instances ithreads odir
	for m in "${_mapper_addreadgroup[@]}"; do
		declare -n _bams_addreadgroup=$m
		((instances+=${#_bams_addreadgroup[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_addreadgroup[@]}"; do
		declare -n _bams_addreadgroup=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_addreadgroup[@]}"; do
			[[ ${nidx[$i]} ]] && rgprefix=NORMAL_ || rgprefix=TUMOR_
			[[ ! $tidx ]] && rgprefix=''

			tomerge=()
			o="$(basename "${_bams_addreadgroup[$i]}")"
			o="${o%.*}"
			while read -r slice; do
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					picard
						-Xmx${jmem}m
						-XX:ParallelGCThreads=$jgct
						-XX:ConcGCThreads=$jcgct
						-Djava.io.tmpdir="$tmpdir"
						AddOrReplaceReadGroups
						I="$slice"
						O="$slice.rg"
						RGID=${rgprefix}1
						RGLB=${rgprefix}lib
						RGPL=illumina
						RGPU=${rgprefix}unit
						RGSM="$o"
						VALIDATION_STRINGENCY=SILENT
						VERBOSITY=WARNING
				CMD
					mv "$slice.rg" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_addreadgroup[${_bams_addreadgroup[$i]}]}"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$odir/$o.rg.bam"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD
			_bamslices_addreadgroup["$odir/$o.rg.bam"]="${_bamslices_addreadgroup[${_bams_addreadgroup[$i]}]}"
			_bams_addreadgroup[$i]="$odir/$o.rg.bam"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t $minstances -a cmd1 && \
			commander::runcmd -v -b -t $instances -a cmd2
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::splitncigar() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_splitncigar _bamslices_splitncigar
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_splitncigar=$OPTARG;;
			c) ((mandatory++)); _bamslices_splitncigar=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	commander::print "splitting N-cigar alignments"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice instances ithreads odir tdir dinstances djmem djgct djcgct
	for m in "${_mapper_splitncigar[@]}"; do
		declare -n _bams_splitncigar=$m
		((instances+=${#_bams_splitncigar[@]}))
	done
	read -r dinstances ithreads djmem djgct djcgct < <(configure::jvm -i $instances -t 1 -T $threads)
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_splitncigar[@]}"; do
		declare -n _bams_splitncigar=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_splitncigar[@]}"; do
			tomerge=()
			o="$(basename "${_bams_splitncigar[$i]}")"
			o="${o%.*}"

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				ln -sfn "$genome" "$tdir/$o.fa"
			CMD
				samtools faidx "$tdir/$o.fa"
			CMD
				rm -f "$tdir/$o.dict"
			CMD
				picard
					-Xmx${djmem}m
					-XX:ParallelGCThreads=$djgct
					-XX:ConcGCThreads=$djcgct
					-Djava.io.tmpdir="$tmpdir"
					CreateSequenceDictionary
					R="$tdir/$o.fa"
					VERBOSITY=WARNING
			CMD

			# v3.X ReassignOneMappingQuality: e.g. misused 255 as unique flag to 60
			# -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
			# v4.X does it automatically for 255, but better do it via STAR parameter outSAMmapqUnique
			# ALLOW_N_CIGAR_READS removed
			# maybe test later: --read-filter NonChimericOriginalAlignmentReadFilter
			# 					--read-filter MappingQualityReadFilter --minimum-mapping-quality 0
			# only for FR-PE data else produces empty file!	--read-filter MateDifferentStrandReadFilter
			while read -r slice; do
				commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					gatk
						--java-options '
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
						'
						SplitNCigarReads
						-I "$slice"
						-O "$slice.nsplit"
						-R "$tdir/$o.fa"
						-verbosity ERROR
						--tmp-dir $tmpdir
				CMD
					mv "$slice.nsplit" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_splitncigar[${_bams_splitncigar[$i]}]}"

			o="$odir/$o.nsplit.bam"
			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_splitncigar["$o"]="${_bamslices_splitncigar[${_bams_splitncigar[$i]}]}"
			_bams_splitncigar[$i]="$o"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	} || {
		{	commander::runcmd -v -b -t $dinstances -a cmd1 && \
            commander::runcmd -v -b -t $minstances -a cmd2 && \
			commander::runcmd -v -b -t $instances -a cmd3
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::leftalign() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_leftalign _bamslices_leftalign
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_leftalign=$OPTARG;;
			c) ((mandatory++)); _bamslices_leftalign=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	commander::print "leftaligning alignments"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir tdir instances ithreads dinstances djmem djgct djcgct
	for m in "${_mapper_leftalign[@]}"; do
		declare -n _bams_leftalign=$m
		((instances+=${#_bams_leftalign[@]}))
	done
   	read -r dinstances ithreads djmem djgct djcgct < <(configure::jvm -i $instances -t 1 -T $threads)
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_leftalign[@]}"; do
		declare -n _bams_leftalign=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_leftalign[@]}"; do
			tomerge=()
			o="$(basename "${_bams_leftalign[$i]}")"
			o="${o%.*}"

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				ln -sfn "$genome" "$tdir/$o.fa"
			CMD
				samtools faidx "$tdir/$o.fa"
			CMD
				rm -f "$tdir/$o.dict"
			CMD
				picard
					-Xmx${djmem}m
					-XX:ParallelGCThreads=$djgct
					-XX:ConcGCThreads=$djcgct
					-Djava.io.tmpdir="$tmpdir"
					CreateSequenceDictionary
					R="$tdir/$o.fa"
					VERBOSITY=WARNING
			CMD

			while read -r slice; do
				commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						LeftAlignIndels
						-I "$slice"
						-O "$slice.leftaln"
						-R "$tdir/$o.fa"
						-verbosity ERROR
						--tmp-dir $tmpdir
				CMD
					mv "$slice.leftaln" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_leftalign[${_bams_leftalign[$i]}]}"

			o="$odir/$o.leftaln.bam"
			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_leftalign["$o"]="${_bamslices_leftalign[${_bams_leftalign[$i]}]}"
			_bams_leftalign[$i]="$o"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
        commander::printcmd -a cmd3
	} || {
		{	commander::runcmd -v -b -t $dinstances -a cmd1 && \
            commander::runcmd -v -b -t $minstances -a cmd2 && \
			commander::runcmd -v -b -t $instances -a cmd3
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::bqsr() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_bqsr _bamslices_bqsr
	declare -A nidx tidx
	while getopts 'S:s:t:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((mandatory++)); _mapper_bqsr=$OPTARG;;
			c) ((mandatory++)); _bamslices_bqsr=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1
	if [[ ! $dbsnp ]]; then
		dbsnp="$tmpdir/$(basename "$genome").vcf"
		echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$dbsnp"
		bgzip -f -@ $threads < "$dbsnp" > "$dbsnp.gz"
		tabix -f -p vcf "$dbsnp.gz"
		dbsnp="$dbsnp.gz"
	fi
	commander::print "base quality score recalibration"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir instances ithreads dinstances djmem djgct djcgct
	for m in "${_mapper_bqsr[@]}"; do
		declare -n _bams_bqsr=$m
		((instances+=${#_bams_bqsr[@]}))
	done
    read -r dinstances ithreads djmem djgct djcgct < <(configure::jvm -i $instances -t 1 -T $threads)
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4
	for m in "${_mapper_bqsr[@]}"; do
		declare -n _bams_bqsr=$m
		odir="$outdir/$m"
        tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_bqsr[@]}"; do
			tomerge=()
			o="$(basename "${_bams_bqsr[$i]}")"
			o="${o%.*}"

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				ln -sfn "$genome" "$tdir/$o.fa"
			CMD
				samtools faidx "$tdir/$o.fa"
			CMD
				rm -f "$tdir/$o.dict"
			CMD
				picard
					-Xmx${djmem}m
					-XX:ParallelGCThreads=$djgct
					-XX:ConcGCThreads=$djcgct
					-Djava.io.tmpdir="$tmpdir"
					CreateSequenceDictionary
					R="$tdir/$o.fa"
					VERBOSITY=WARNING
			CMD

			while read -r slice; do
				# https://gatkforums.broadinstitute.org/gatk/discussion/7131/is-indel-realignment-removed-from-gatk4
				# https://software.broadinstitute.org/gatk/blog?id=7847
				#
				# known polymorphic sites used to exclude regions around
				# https://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x
				# gatk known-sites bundle: https://software.broadinstitute.org/gatk/documentation/article.php?id=1213
				#
				# BaseRecalibratorSpark with --spark-master local[$mthreads] is BETA!!
				# --bqsr-baq-gap-open-penalty (Phred Scaled) Default value is 40. 30 is perhaps better for whole genome call sets
				# ApplyBQSRSpark fails as of v4.1.2.0
				# -> i.e. also true for BQSRPipelineSpark which does both steps in one
				# GatherBQSRReports - Gathers scattered BQSR recalibration reports into a single file
				# 
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						BaseRecalibrator
						-R "$tdir/$o.fa"
						--known-sites "$dbsnp"
						-I "$slice"
						-O "$slice.bqsreport"
						-verbosity ERROR
						--tmp-dir $tmpdir
				CMD

				commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					rm -rf "$slice.bqsr.parts"
				CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						ApplyBQSR
						-bqsr "$slice.bqsreport"
						-I "$slice"
						-O "$slice.bqsr"
						-verbosity ERROR
						--tmp-dir $tmpdir
				CMD
					mv "$slice.bqsr" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_bqsr[${_bams_bqsr[$i]}]}"

			o="$odir/$o.bqsr.bam"
			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_bqsr["$o"]="${_bamslices_bqsr[${_bams_bqsr[$i]}]}"
			_bams_bqsr[$i]="$o"
		done
	done

	$skip && {
        commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	} || {
		{	commander::runcmd -v -b -t $dinstances -a cmd1 && \
            commander::runcmd -v -b -t $minstances -a cmd2 && \
			commander::runcmd -v -b -t $minstances -a cmd3 && \
			commander::runcmd -v -b -t $instances -a cmd4
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
