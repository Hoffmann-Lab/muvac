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
		{	commander::runcmd -v -b -t 1 -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd1
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

callvariants::haplotypecaller() {
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
	declare -n _mapper_haplotypecaller _bamslices_haplotypecaller
	declare -A nidx tidx
	while getopts 'S:s:t:g:d:m:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((mandatory++)); _mapper_haplotypecaller=$OPTARG;;
			c) ((mandatory++)); _bamslices_haplotypecaller=$OPTARG;;
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

	commander::print "calling variants haplotypecaller"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o e slice odir instances ithreads
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		((instances+=${#_bams_haplotypecaller[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_haplotypecaller[@]}"; do
			tomerge=()
            
			while read -r slice; do
				# all alt Phredscaled Likelihoods ordering:
				# reference homozygous (0/0)
				# ref and alt 1 heterzygous (0/1)
				# alt 1 homozygous (1/1)
				# ref and alt 2 heterzygous (0/2)
				# alt 1 and alt 2 heterozygous (1/2)
				# alt 2 homozygous (2/2)
				# gatk bug as of v4.1.2.0 --max-reads-per-alignment-start 0 not a valid option
				# HaplotypeCallerSpark with --spark-master local[$mthreads] is BETA and differs in results!!
				# Spark does not yet support -D "$dbsnp"
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						HaplotypeCaller
						-I "$slice"
						-O "$slice.vcf"
						-R "$genome"
						-D "$dbsnp"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--min-base-quality-score 20
						--native-pair-hmm-threads $mthreads
						--max-alternate-alleles 3
						--all-site-pls true
						-verbosity INFO
						--tmp-dir $tmpdir
				CMD

				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$slice.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_haplotypecaller[${_bams_haplotypecaller[$i]}]}"

			o="$odir/$(basename "${_bams_haplotypecaller[$i]}")"
			o="${o%.*}"

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$o.vcf" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T $tmpdir -m ${memory}M -o "$o.$e"
				CMD

				commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
	} || {
		{	commander::runcmd -v -b -t $minstances -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			commander::runcmd -v -b -t $threads -a cmd4 && \
			commander::runcmd -v -b -t $minstances -a cmd5 && \
			commander::runcmd -v -b -t $instances -a cmd6
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

callvariants::mutect() {
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
			-1 <normalidx> | array of (requieres -2)
			-2 <tumoridx>  | array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_mutect _bamslices_mutect _nidx_mutect _tidx_mutect
	declare -A nidx tidx
	while getopts 'S:s:t:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_mutect=$OPTARG;;
			c) ((mandatory++)); _bamslices_mutect=$OPTARG;;
			1) ((mandatory++)); _nidx_mutect=$OPTARG;;
			2) ((mandatory++)); _tidx_mutect=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage && return 1

	commander::print "calling variants mutect"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir ithreads instances=$((${#_mapper_mutect[@]}*${#_tidx_mutect[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_mutect[@]}"; do
		declare -n _bams_mutect=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_tidx_mutect[@]}"; do
			tomerge=()
            
			while read -r nslice slice; do
				# normal name defined for RGSM sam header entry by alignment::addreadgroup 
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						Mutect2
						-I "$nslice"
						-I "$slice"
						-normal NORMAL
						-tumor TUMOR
						-O "$slice.unsorted.vcf"
						-R "$genome"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--min-base-quality-score 20
						--native-pair-hmm-threads $mthreads
						-verbosity INFO
						--tmp-dir $tmpdir
				CMD

				# commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				# 	vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				# CMD

				# commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				# 	bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				# CMD
				# 	vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				# CMD

				# commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				# 	vcffixup "$slice.fixed.nomulti.vcf"
				# CMD
				# 	vt normalize -q -n -r "$genome" -
				# CMD
				# 	vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
				# CMD

				tomerge+=("$slice")
			done < <(paste "${_bamslices_mutect[${_bams_mutect[${_nidx_mutect[$i]}]}]}" "${_bamslices_mutect[${_bams_mutect[${_tidx_mutect[$i]}]}]}")

			o="$odir/$(basename "${_bams_mutect[${_tidx_mutect[$i]}]}")"
			o="${o%.*}"

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$o.vcf" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T $tmpdir -m ${memory}M -o "$o.$e"
				CMD

				commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD

				break
			done

		done
	done

	$skip && {
		commander::printcmd -a cmd1
	# 	commander::printcmd -a cmd2
	# 	commander::printcmd -a cmd3
	# 	commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
	} || {
		{	commander::runcmd -v -b -t $minstances -a cmd1 && \
	# 		commander::runcmd -v -b -t $threads -a cmd2 && \
	# 		commander::runcmd -v -b -t $threads -a cmd3 && \
	# 		commander::runcmd -v -b -t $threads -a cmd4 && \
			commander::runcmd -v -b -t $minstances -a cmd5 && \
			commander::runcmd -v -b -t $instances -a cmd6
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
