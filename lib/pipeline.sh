#! /usr/bin/env bash
# (c) Konstantin Riege

pipeline::index(){
	unset NA1 NA2
	alignment::segemehl \
		-S ${NOsege:=false} \
		-s true \
		-t $THREADS \
		-g "$GENOME" \
		-x "$GENOME.segemehl.idx" \
		-o "$TMPDIR" \
		-r NA1 \
		-1 NA2
	unset NA1 NA2
	alignment::star \
		-S ${NOstar:=false} \
		-s true \
		-t $THREADS \
		-g "$GENOME" \
		-x "$GENOME.star.idx" \
		-f "$GTF" \
		-o "$TMPDIR" \
		-p "$TMPDIR" \
		-r NA1 \
		-1 NA2
	unset NA1 NA2
	alignment::bwa \
		-S ${NObwa:=false} \
		-s true \
		-t $THREADS \
		-g "$GENOME" \
		-x "$GENOME.bwa.idx/bwa" \
		-o "$TMPDIR" \
		-r NA1 \
		-1 NA2
	genome::mkdict \
		-t $THREADS \
		-i "$GENOME" \
		-p "$TMPDIR"
	return 0
}

pipeline::_slice(){
	alignment::slice \
		-S ${SLICED:-$1} \
		-s ${SKIPslice:-$2} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" || return 1
	! $1 && ! $2 && SLICED=true
	! $1 && ${SKIPslice:-false} && SLICED=true

	return 0
}

pipeline::_preprocess(){
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/raw" -1 FASTQ1 -2 FASTQ2
		preprocess::fastqc \
			-S ${NOqual:=false} \
			-s ${SKIPfqual:=false} \
			-t $THREADS \
			-o "$OUTDIR/qualities/raw" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2

		${NOtrim:=false} || {
			preprocess::trimmomatic \
				-S ${NOtrim:=false} \
				-s ${SKIPtrim:=false} \
				-t $THREADS \
				-o "$OUTDIR/trimmed" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/trimmed" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPfqual:=false} \
				-t $THREADS \
				-o "$OUTDIR/qualities/trimmed" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		${NOclip:=false} || {
			preprocess::rmpolynt \
				-S ${NOclip:=false} \
				-s ${SKIPclip:=false} \
				-t $THREADS \
				-o "$OUTDIR/polyntclipped" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/polyntclipped" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPfqual:=false} \
				-t $THREADS \
				-o "$OUTDIR/qualities/polyntclipped" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		if [[ $ADAPTER1 ]]; then
			${NOclip:=false} || {
				preprocess::cutadapt \
					-S ${NOclip:=false} \
					-s ${SKIPclip:=false} \
					-a ADAPTER1 \
					-A ADAPTER2 \
					-t $THREADS \
					-o "$OUTDIR/adapterclipped" \
					-1 FASTQ1 \
					-2 FASTQ2
				preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/adapterclipped" -1 FASTQ1 -2 FASTQ2
				preprocess::fastqc \
					-S ${NOqual:=false} \
					-s ${SKIPfqual:=false} \
					-t $THREADS \
					-o "$OUTDIR/qualities/adapterclipped" \
					-p "$TMPDIR" \
					-1 FASTQ1 \
					-2 FASTQ2
			}
		fi

		preprocess::rcorrector \
			-S ${NOcor:=true} \
			-s ${SKIPcor:=false} \
			-t $THREADS \
			-o "$OUTDIR/corrected" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2

		${NOrrm:=true} || {
			preprocess::sortmerna \
				-S ${NOrrm:=true} \
				-s ${SKIPrrm:=false} \
				-t $THREADS \
				-o "$OUTDIR/rrnafiltered" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/rrnafiltered" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPfqual:=false} \
				-t $THREADS \
				-o "$OUTDIR/qualities/rrnafiltered" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::qcstats \
			-S ${NOstats:=false} \
			-s ${SKIPstats:=false} \
			-t $THREADS \
			-i qualdirs \
			-o "$OUTDIR/stats" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2
	fi

	return 0
}


pipeline::_mapping(){
	if [[ ! $MAPPED ]]; then
		alignment::segemehl \
			-S ${NOsege:=false} \
			-s ${SKIPsege:=false} \
			-5 ${SKIPmd5:=false} \
			-1 FASTQ1 \
			-2 FASTQ2 \
			-o "$OUTDIR/mapped" \
			-t $THREADS \
			-a $((100-DISTANCE)) \
			-i ${INSERTSIZE:=200000} \
			-n ${NOsplitreads:=true} \
			-g "$GENOME" \
			-x "$GENOME.segemehl.idx" \
			-r mapper

		alignment::star \
			-S ${NOstar:=false} \
			-s ${SKIPstar:=false} \
			-5 ${SKIPmd5:=false} \
			-1 FASTQ1 \
			-2 FASTQ2 \
			-o "$OUTDIR/mapped" \
			-p "$TMPDIR" \
			-t $THREADS \
			-a $((100-DISTANCE)) \
			-i ${INSERTSIZE:=200000} \
			-n ${NOsplitreads:=true} \
			-g "$GENOME" \
			-f "$GTF" \
			-x "$GENOME.star.idx" \
			-r mapper

		! ${NOsplitreads:=true} || alignment::bwa \
			-S ${NObwa:=false} \
			-s ${SKIPbwa:=false} \
			-5 ${SKIPmd5:=false} \
			-1 FASTQ1 \
			-2 FASTQ2 \
			-o "$OUTDIR/mapped" \
			-t $THREADS \
			-a $((100-DISTANCE)) \
			-f true \
			-g "$GENOME" \
			-x "$GENOME.bwa.idx/bwa" \
			-r mapper
	else
		declare -g -a ${MAPNAME:=custom}
		declare -n _MAPNAME_muvac=$MAPNAME
		_MAPNAME_muvac=("${MAPPED[@]}")
		mapper+=($MAPNAME)
	fi

	[[ ${#mapper[@]} -eq 0 ]] && return 0

	alignment::add4stats -r mapper
	alignment::bamqc \
		-S ${NOqual:=false} \
		-s ${SKIPmqual:=false} \
		-t $THREADS \
		-r mapper

	alignment::postprocess \
		-S ${NOuniq:=false} \
		-s ${SKIPuniq:=false} \
		-j uniqify \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper
	! ${NOuniq:=false} && ${NOsort:=false} && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::postprocess \
		-S ${NOsort:=false} \
		-s ${SKIPsort:=false} \
		-j sort \
		-t "$THREADS" \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper
	(${NOuniq:=false} && ! ${NOsort:=false}) || (! ${NOuniq:=false} && ! ${NOsort:=false}) && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	return 0
}

pipeline::germline() {
	declare -a mapper
	declare -A slicesinfo

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${NOdict:=false} \
		-s ${SKIPdict:=false} \
		-5 ${SKIPmd5:=false} \
		-i "$GENOME" \
		-p "$TMPDIR" \
		-t $THREADS

	pipeline::_slice ${NOaddrg:=false} ${SKIPaddrg:=false}
	alignment::addreadgroup \
		-S ${NOaddrg:=false} \
		-s ${SKIPaddrg:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-n ${RGPREFIX:='SAMPLE'} \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	pipeline::_slice ${NOrmd:=false} ${SKIPrmd:=false}
	${NOrmd:=false} || {
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${NOcmo:=false} ${SKIPcmo:=false}
	${NOcmo:=false} || {
		alignment::clipmateoverlaps \
			-S ${NOcmo:=false} \
			-s ${SKIPcmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::qcstats \
		-S ${NOstats:=false} \
		-s ${SKIPstats:=false} \
		-r mapper \
		-t $THREADS \
		-o "$OUTDIR/stats"

	pipeline::_slice ${NOreo:=true} ${SKIPreo:=false}
	alignment::reorder \
		-S ${NOreo:=true} \
		-s ${SKIPreo:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	${NOsplitreads:=true} || {
		pipeline::_slice ${NOnsplit:=true} ${SKIPnsplit:=false}
		alignment::splitncigar \
			-S ${NOnsplit:=false} \
			-s ${SKIPnsplit:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-g "$GENOME" \
			-r mapper \
			-c slicesinfo \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
	}

	pipeline::_slice ${NOlaln:=false} ${SKIPlaln:=false}
	alignment::leftalign \
		-S ${NOlaln:=false} \
		-s ${SKIPlaln:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	! ${NOdbsnp:-false} && [[ $DBSNP ]] && {
		variants::vcfzip -t $THREADS -i "$DBSNP"
	}

	pipeline::_slice ${NObqsr:=false} ${SKIPbqsr:=false}
	alignment::bqsr \
		-S ${NObqsr:=false} \
		-s ${SKIPbqsr:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g $GENOME \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	alignment::postprocess \
		-S ${NOidx:=false} \
		-s ${SKIPidx:=false} \
		-j index \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper

	${NOpon:=true} || {
		pipeline::_slice ${NOpon:=true} ${SKIPpon:=false}
		variants::panelofnormals \
			-S ${NOpon:=true} \
			-s ${SKIPpon:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-p "$TMPDIR" \
			-o "$OUTDIR/pon"

		variants::makepondb \
			-S ${NOpondb:=false} \
			-s ${SKIPpondb:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-M $MAXMEMORY \
			-r mapper \
			-p "$TMPDIR" \
			-o "$OUTDIR/pon"

		return 0
	}

	pipeline::_slice ${NOgatk:=false} ${SKIPgatk:=false}
	variants::haplotypecaller \
		-S ${NOgatk:=false} \
		-s ${SKIPgatk:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-e ${NOsplitreads:=true} \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::bcftools \
		-S ${NObt:=false} \
		-s ${SKIPbt:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::freebayes \
		-S ${NOfb:=false} \
		-s ${SKIPfb:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::varscan \
		-S ${NOvs:=false} \
		-s ${SKIPvs:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::vardict \
		-S ${NOvd:=false} \
		-s ${SKIPvd:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::platypus \
		-S ${NOpp:=false} \
		-s ${SKIPpp:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	return 0
}

pipeline::somatic() {
	declare -a mapper
	declare -A slicesinfo

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${NOdict:=false} \
		-s ${SKIPdict:=false} \
		-5 ${SKIPmd5:=false} \
		-i "$GENOME" \
		-p "$TMPDIR" \
		-t $THREADS

	pipeline::_slice ${NOaddrg:=false} ${SKIPaddrg:=false}
	alignment::addreadgroup \
		-S ${NOaddrg:=false} \
		-s ${SKIPaddrg:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	pipeline::_slice ${NOrmd:=false} ${SKIPrmd:=false}
	${NOrmd:=false} || {
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${NOcmo:=false} ${SKIPcmo:=false}
	${NOcmo:=false} || {
		alignment::clipmateoverlaps \
			-S ${NOcmo:=false} \
			-s ${SKIPcmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPmqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::qcstats \
		-S ${NOstats:=false} \
		-s ${SKIPstats:=false} \
		-r mapper \
		-t $THREADS \
		-o "$OUTDIR/stats"

	pipeline::_slice ${NOreo:=true} ${SKIPreo:=false}
	alignment::reorder \
		-S ${NOreo:=true} \
		-s ${SKIPreo:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	${NOsplitreads:=true} || {
		pipeline::_slice ${NOnsplit:=true} ${SKIPnsplit:=false}
		alignment::splitncigar \
			-S ${NOnsplit:=false} \
			-s ${SKIPnsplit:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-g "$GENOME" \
			-r mapper \
			-c slicesinfo \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
	}

	pipeline::_slice ${NOlaln:=false} ${SKIPlaln:=false}
	alignment::leftalign \
		-S ${NOlaln:=false} \
		-s ${SKIPlaln:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	! ${NOdbsnp:-false} && [[ $DBSNP ]] && {
		variants::vcfzip -t $THREADS -i "$DBSNP"
	}

	pipeline::_slice ${NObqsr:=false} ${SKIPbqsr:=false}
	alignment::bqsr \
		-S ${NObqsr:=false} \
		-s ${SKIPbqsr:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped"

	alignment::postprocess \
		-S ${NOidx:=false} \
		-s ${SKIPidx:=false} \
		-j index \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper

	pipeline::_slice ${NOgatk:=false} ${SKIPgatk:=false}
	variants::mutect \
		-S ${NOgatk:=false} \
		-s ${SKIPgatk:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-n "$(${NOpon:-false} || echo "$PONDB")" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-c slicesinfo \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::bcftools \
		-S ${NObt:=false} \
		-s ${SKIPbt:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::freebayes \
		-S ${NOfb:=false} \
		-s ${SKIPfb:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::varscan \
		-S ${NOvs:=false} \
		-s ${SKIPvs:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::vardict \
		-S ${NOvd:=false} \
		-s ${SKIPvd:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	variants::platypus \
		-S ${NOpp:=false} \
		-s ${SKIPpp:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$(${NOdbsnp:-false} || echo "$DBSNP")" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-p "$TMPDIR" \
		-o "$OUTDIR/variants"

	return 0
}
