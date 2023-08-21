#! /usr/bin/env bash
# (c) Konstantin Riege

function pipeline::index(){
	genome::mkdict \
		-F \
		-t $THREADS \
		-i "$GENOME"

	unset NA1 NA2
	alignment::segemehl \
		-S ${NOsege:=false} \
		-s true \
		-t $THREADS \
		-g "$GENOME" \
		-x "$GENOME.segemehl.idx" \
		-o "$TMPDIR" \
		-F \
		-r NA1 \
		-1 NA2
	unset NA1 NA2
	alignment::star \
		-S ${NOstar:=false} \
		-s true \
		-t $THREADS \
		-g "$GENOME" \
		-x "$GENOME.star.idx" \
		-o "$TMPDIR" \
		-F \
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
		-F \
		-r NA1 \
		-1 NA2

	return 0
}

function pipeline::_slice(){
	alignment::slice \
		-S ${SLICED:-$1} \
		-s ${SKIPslice:-$2} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-r mapper \
		-c slicesinfo
	! $1 && ! $2 && SLICED=true
	! $1 && ${SKIPslice:-false} && SLICED=true

	return 0
}

function pipeline::_preprocess(){
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		local params=""
		[[ $ADAPTER1 ]] || params="-a ADAPTER1 -A ADAPTER2"
		preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/raw" -1 FASTQ1 -2 FASTQ2
		preprocess::fastqc \
			-S ${NOqual:=false} \
			-s ${SKIPfqual:=false} \
			-t $THREADS \
			-M $MAXMEMORY \
			-o "$OUTDIR/qualities/raw" \
			-1 FASTQ1 \
			-2 FASTQ2 \
			$params

		${NOtrim:=false} || {
			preprocess::trimmomatic \
				-S ${NOtrim:=false} \
				-s ${SKIPtrim:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/trimmed" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/trimmed" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPtrim:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/trimmed" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		${NOpclip:=true} || {
			preprocess::rmpolynt \
				-S ${NOpclip:=true} \
				-s ${SKIPpclip:=false} \
				-t $THREADS \
				-o "$OUTDIR/polyntclipped" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/polyntclipped" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPpclip:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/polyntclipped" \
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
					-s ${SKIPclip:=false} \
					-t $THREADS \
					-M $MAXMEMORY \
					-o "$OUTDIR/qualities/adapterclipped" \
					-1 FASTQ1 \
					-2 FASTQ2
			}
		fi

		preprocess::rcorrector \
			-S ${NOcor:=true} \
			-s ${SKIPcor:=false} \
			-t $THREADS \
			-o "$OUTDIR/corrected" \
			-1 FASTQ1 \
			-2 FASTQ2

		${NOrrm:=true} || {
			preprocess::sortmerna \
				-S ${NOrrm:=true} \
				-s ${SKIPrrm:=false} \
				-t $THREADS \
				-o "$OUTDIR/rrnafiltered" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/rrnafiltered" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPrrm:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/rrnafiltered" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::qcstats \
			-S ${NOstats:=false} \
			-s ${SKIPstats:=false} \
			-t $THREADS \
			-i qualdirs \
			-o "$OUTDIR/stats" \
			-1 FASTQ1 \
			-2 FASTQ2
	fi

	return 0
}


function pipeline::_mapping(){
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
			-t $THREADS \
			-a $((100-DISTANCE)) \
			-i ${INSERTSIZE:=200000} \
			-n ${NOsplitreads:=true} \
			-g "$GENOME" \
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
		-o "$OUTDIR/mapped" \
		-r mapper
	! ${NOuniq:=false} && ${NOsort:=false} && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPuniq:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::postprocess \
		-S ${NOsort:=false} \
		-s ${SKIPsort:=false} \
		-j sort \
		-t "$THREADS" \
		-o "$OUTDIR/mapped" \
		-r mapper
	alignment::postprocess \
		-S ${NOsort:=false} \
		-s ${SKIPsort:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

	(${NOuniq:=false} && ! ${NOsort:=false}) || (! ${NOuniq:=false} && ! ${NOsort:=false}) && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPsort:=false} \
			-t $THREADS \
			-r mapper
	}

	return 0
}

function pipeline::germline(){
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOaddrg:=false} \
		-s ${SKIPaddrg:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

	pipeline::_slice ${NOrmd:=false} ${SKIPrmd:=false}
	${NOrmd:=false} || {
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-3 FASTQ3 \
			-c slicesinfo \
			-x "$REGEX" \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPrmd:=false} \
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
		alignment::postprocess \
			-S ${NOcmo:=false} \
			-s ${SKIPcmo:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPcmo:=false} \
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOreo:=true} \
		-s ${SKIPreo:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

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
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${NOnsplit:=false} \
			-s ${SKIPnsplit:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOlaln:=false} \
		-s ${SKIPlaln:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

	pipeline::_slice ${NObqsr:=false} ${SKIPbqsr:=false}
	alignment::bqsr \
		-S ${NObqsr:=false} \
		-s ${SKIPbqsr:=false} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-g $GENOME \
		-d "$DBSNP" \
		-r mapper \
		-c slicesinfo \
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NObqsr:=false} \
		-s ${SKIPbqsr:=false} \
		-j index \
		-t $THREADS \
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
			-o "$OUTDIR/pon"

		variants::makepondb \
			-S ${NOpondb:=false} \
			-s ${SKIPpondb:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-M $MAXMEMORY \
			-r mapper \
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
		-d "$DBSNP" \
		-r mapper \
		-e ${NOsplitreads:=true} \
		-c slicesinfo \
		-o "$OUTDIR/variants"

	variants::bcftools \
		-S ${NObt:=false} \
		-s ${SKIPbt:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-o "$OUTDIR/variants"

	variants::freebayes \
		-S ${NOfb:=false} \
		-s ${SKIPfb:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-o "$OUTDIR/variants"

	variants::varscan \
		-S ${NOvs:=false} \
		-s ${SKIPvs:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-o "$OUTDIR/variants"

	variants::vardict \
		-S ${NOvd:=false} \
		-s ${SKIPvd:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-o "$OUTDIR/variants"

	variants::platypus \
		-S ${NOpp:=false} \
		-s ${SKIPpp:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-o "$OUTDIR/variants"

	return 0
}

function pipeline::somatic(){
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOaddrg:=false} \
		-s ${SKIPaddrg:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

	pipeline::_slice ${NOrmd:=false} ${SKIPrmd:=false}
	${NOrmd:=false} || {
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-r FASTQ3 \
			-c slicesinfo \
			-x "$REGEX" \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPrmd:=false} \
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
		alignment::postprocess \
			-S ${NOcmo:=false} \
			-s ${SKIPcmo:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${NOqual:=false} \
			-s ${SKIPcmo:=false} \
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOreo:=true} \
		-s ${SKIPreo:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

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
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${NOnsplit:=false} \
			-s ${SKIPnsplit:=false} \
			-j index \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper
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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NOlaln:=false} \
		-s ${SKIPlaln:=false} \
		-j index \
		-t $THREADS \
		-o "$OUTDIR/mapped" \
		-r mapper

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
		-o "$OUTDIR/mapped"
	alignment::postprocess \
		-S ${NObqsr:=false} \
		-s ${SKIPbqsr:=false} \
		-j index \
		-t $THREADS \
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
		-n "$PONDB" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-c slicesinfo \
		-o "$OUTDIR/variants"

	variants::bcftools \
		-S ${NObt:=false} \
		-s ${SKIPbt:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-o "$OUTDIR/variants"

	variants::freebayes \
		-S ${NOfb:=false} \
		-s ${SKIPfb:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-o "$OUTDIR/variants"

	variants::varscan \
		-S ${NOvs:=false} \
		-s ${SKIPvs:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-o "$OUTDIR/variants"

	variants::vardict \
		-S ${NOvd:=false} \
		-s ${SKIPvd:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-o "$OUTDIR/variants"

	variants::platypus \
		-S ${NOpp:=false} \
		-s ${SKIPpp:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-d "$DBSNP" \
		-r mapper \
		-1 NIDX \
		-2 TIDX \
		-o "$OUTDIR/variants"

	return 0
}
