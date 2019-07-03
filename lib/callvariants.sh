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

	commander::print "calling variants haplotypecaller"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir instances ithreads dinstances djmem djgct djcgct
	for m in "${_mapper_bqsr[@]}"; do
		declare -n _bams_bqsr=$m
		((instances+=${#_bams_bqsr[@]}))
	done
	read -r dinstances ithreads djmem djgct djcgct < <(configure::jvm -i $instances -t 1 -T $threads)
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
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
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
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
						-R "$tdir/$o.fa"
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

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$slice.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_bqsr[${_bams_bqsr[$i]}]}"

			o="$odir/$o"

			commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				bcftools concat -o "$o.vcf" "${tomerge[@]/%/.vcf}"
			CMD
				bcftools concat --threads $ithreads -O z -o "$o.vcf.gz" "${tomerge[@]/%/.vcf}"
			CMD
				tabix -f -p vcf "$o.vcf.gz"
			CMD
			commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				bcftools concat -o "$o.fixed.vcf" "${tomerge[@]/%/.fixed.vcf}"
			CMD
				bcftools concat --threads $ithreads -O z -o "$o.fixed.vcf.gz" "${tomerge[@]/%/.fixed.vcf}"
			CMD
				tabix -f -p vcf "$o.fixed.vcf.gz"
			CMD
			commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				bcftools concat -o "$o.fixed.nomulti.vcf" "${tomerge[@]/%/.fixed.nomulti.vcf}"
			CMD
				bcftools concat --threads $ithreads -O z -o "$o.fixed.nomulti.vcf.gz" "${tomerge[@]/%/.fixed.nomulti.vcf}"
			CMD
				tabix -f -p vcf "$o.fixed.nomulti.vcf.gz"
			CMD
			commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				bcftools concat -o "o.fixed.nomulti.normed.vcf" "${tomerge[@]/%/.fixed.nomulti.normed.vcf}"
			CMD
				bcftools concat --threads $ithreads -O z -o "$o.fixed.nomulti.normed.vcf.gz" "${tomerge[@]/%/.fixed.nomulti.normed.vcf}"
			CMD
				tabix -f -p vcf "$o.fixed.nomulti.normed.vcf.gz"
			CMD
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
		{	commander::runcmd -v -b -t $dinstances -a cmd1 && \
            commander::runcmd -v -b -t $minstances -a cmd2 && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			commander::runcmd -v -b -t $threads -a cmd4 && \
			commander::runcmd -v -b -t $threads -a cmd5 && \
			commander::runcmd -v -b -t $instances -a cmd6
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
