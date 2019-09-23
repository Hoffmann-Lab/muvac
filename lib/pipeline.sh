#! /usr/bin/env bash
# (c) Konstantin Riege

pipeline::_preprocess(){
	source $INSDIR/conda/bin/activate py2

	${SKIPmd5:=false} || {
		[[ ! -s $GENOME.md5.sh ]] && cp $INSDIR/latest/bashbone/lib/md5.sh $GENOME.md5.sh
		source $GENOME.md5.sh
		thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
		[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	}
	
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		{	qualdirs+=("$OUTDIR/qualities/raw") && \
			preprocess::fastqc \
				-S ${NOqual:=false} \
				-s ${SKIPqual:=false} \
				-t $THREADS \
				-o $OUTDIR/qualities/raw \
				-1 FASTQ1 \
				-2 FASTQ2
		} || return 1
		if [[ $ADAPTER ]]; then
			${NOclip:=false} || {
				{	qualdirs+=("$OUTDIR/qualities/clipped") && \
					preprocess::cutadapt \
						-S ${NOclip:=false} \
						-s ${SKIPclip:=false} \
						-a ADAPTER \
						-t $THREADS \
						-o $OUTDIR/clipped \
						-1 FASTQ1 \
						-2 FASTQ2 && \
					preprocess::fastqc \
						-S ${NOqual:=false} \
						-s ${SKIPqual:=false} \
						-t $THREADS \
						-o $OUTDIR/qualities/clipped \
						-1 FASTQ1 \
						-2 FASTQ2
				} || return 1
			}
		fi
		${NOtrim:=false} || { 
			{	qualdirs+=("$OUTDIR/qualities/trimmed") && \
				preprocess::trimmomatic \
					-S ${NOtrim:=false} \
					-s ${SKIPtrim:=false} \
					-t $THREADS \
					-m $MEMORY \
					-o $OUTDIR/trimmed \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${NOqual:=false} \
					-s ${SKIPqual:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/trimmed \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${NOcor:=true} || {
			{	# qualdirs+=("$OUTDIR/qualities/corrected")
				preprocess::rcorrector \
					-S ${NOcor:=false} \
					-s ${SKIPcor:=false} \
					-t $THREADS \
					-o $OUTDIR/corrected \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${NOqual:=false} \
					-s ${SKIPqual:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/corrected \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${NOrrm:=true} || {
			{	qualdirs+=("$OUTDIR/qualities/rrnafiltered") && \
				preprocess::sortmerna \
					-S ${NOrrm:=false} \
					-s ${SKIPrrm:=false} \
					-t $THREADS \
					-m $MEMORY \
					-i $INSDIR \
					-o $OUTDIR/rrnafiltered \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${NOqual:=false} \
					-s ${SKIPqual:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/rrnafiltered \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${NOstats:=false} || {
			{	preprocess::qcstats \
					-S ${NOstats:=false} \
					-s ${SKIPstats:=false} \
					-i qualdirs \
					-o $OUTDIR/stats \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${NOsege:=false} || {
			{	alignment::segemehl \
					-S ${NOsege:=false} \
					-s ${SKIPsege:=false} \
					-5 ${SKIPmd5:=false} \
					-1 FASTQ1 \
					-2 FASTQ2 \
					-o $OUTDIR/mapped \
					-t $THREADS \
					-a $((100-DISTANCE)) \
					-i ${INSERTSIZE:=200000} \
					-p ${NOsplitreads:=true} \
					-g $GENOME \
					-x $GENOME.segemehl.idx \
					-r mapper && \
				alignment::add4stats -r mapper
			} || return 1
		}
	else
		custom=("${MAPPED[@]}")
		mapper+=(custom)
	fi

	[[ ${#mapper[@]} -eq 0 ]] && return 0

	{	alignment::postprocess \
			-S ${nouniq:=false} \
			-s ${SKIPuniq:=false} \
			-j uniqify \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::add4stats -r mapper && \
		alignment::postprocess \
			-S ${NOsort:=false} \
			-s ${SKIPsort:=false} \
			-j sort \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper
		# alignment::postprocess \ <- applied by alignment::slice anyways
		# 	-S ${NOidx:=false} \
		# 	-s ${SKIPidx:=false} \
		# 	-j index \
		# 	-t $THREADS \
		# 	-p $TMPDIR \
		# 	-o $OUTDIR/mapped \
		# 	-r mapper
	} || return 1

	return 0
}

pipeline::_slice(){
	alignment::slice \
		-S $sliced \
		-s $(${SKIPslice:=false} && echo true || echo $1) \
		-t $THREADS \
		-m $MEMORY \
		-r mapper \
		-c slicesinfo \
		-p $TMPDIR || return 1

	$1 || sliced=true # i.e. if not skiptool: sliced=true and -S NOslice=true, else just by SKIPslices slicesinfo will be further updated

	return 0
}

pipeline::germline() {
	declare -a mapper
	declare -A slicesinfo

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${NOdict:=false} \
		-s ${SKIPdict:=false} \
		-5 ${SKIPmd5:=false} \
		-i $GENOME \
		-p $TMPDIR \
		-t $THREADS || return 1

	local sliced=false

	{	pipeline::_slice $($sliced || ${SKIPrg:=false} || ${NOrg:=false} && echo true || echo false) && \
		alignment::addreadgroup \
			-S ${NOrg:=false} \
			-s ${SKIPrg:=false} \
			-t $THREADS \
			-m $MEMORY \
			-n ${RGPREFIX:=''} \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		pipeline::_slice $($sliced || ${SKIPrmd:=false} || ${NOrmd:=false} && echo true || echo false) && \
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		alignment::add4stats -r mapper && \
		alignment::bamstats \
			-S ${NOstats:=false} \
			-s ${SKIPstats:=false} \
			-r mapper \
			-t $THREADS \
			-o $OUTDIR/stats
	} || return 1

	${NOsplitreads:=true} || {
		{ 	pipeline::_slice $($sliced || ${SKIPnsplit:=false} || ${NOnsplit:=false} && echo true || echo false) && \
			alignment::splitncigar \
				-S ${NOnsplit:=false} \
				-s ${SKIPnsplit:=false} \
				-t $THREADS \
				-m $MEMORY \
				-g $GENOME \
				-r mapper \
				-c slicesinfo \
				-p $TMPDIR \
				-o $OUTDIR/mapped
		} || return 1
	}

	{	pipeline::_slice $($sliced || ${SKIPreo:=false} || ${NOreo:=false} && echo true || echo false) && \
		alignment::reorder \
			-S ${NOreo:=false} \
			-s ${SKIPreo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		pipeline::_slice $($sliced || ${SKIPlaln:=false} || ${NOlaln:=false} && echo true || echo false) && \
		alignment::leftalign \
			-S ${NOlaln:=false} \
			-s ${SKIPlaln:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped
	} || return 1

	if [[ $DBSNP ]]; then
		callvariants::vcfzip -t $THREADS -i $DBSNP || return 1
	fi

	{	pipeline::_slice $($sliced || ${SKIPbqsr:=false} || ${NObqsr:=false} && echo true || echo false) && \
		alignment::bqsr \
			-S ${NObqsr:=false} \
			-s ${SKIPbqsr:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-d "$(${NOdbsnp:-false} || echo $DBSNP)" \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		alignment::postprocess \
			-S ${NOidx:=false} \
			-s ${SKIPidx:=false} \
			-j index \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		pipeline::_slice $($sliced || ${SKIPhc:=false} || ${NOhc:=false} && echo true || echo false) && \
		callvariants::haplotypecaller \
			-S ${NOhc:=false} \
			-s ${SKIPhc:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-d "$(${NOdbsnp:-false} || echo $DBSNP)" \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/variants
	} || return 1

	return 0

	{	variants::haplotypecaller && \
		variants::mutect && \
		variants::samtools && \
		variants::varscan && \
		variants::lofreq && \
		variants::platypus && \
		variants::vardict && \
		variants::freebayes && \
		variants::condensevcf && \
		variants::mergevcf
	} || return 1

	return 0
}

pipeline::somatic() {
	declare -a mapper
	declare -A slicesinfo

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${NOdict:=false} \
		-s ${SKIPdict:=false} \
		-5 ${SKIPmd5:=false} \
		-i $GENOME \
		-p $TMPDIR \
		-t $THREADS || return 1

	local sliced=false

	{	pipeline::_slice $($sliced || ${SKIPrg:=false} || ${NOrg:=false} && echo true || echo false) && \
		alignment::addreadgroup \
			-S ${NOrg:=false} \
			-s ${SKIPrg:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-1 NIDX \
			-2 TIDX \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		pipeline::_slice $($sliced || ${SKIPrmd:=false} || ${NOrmd:=false} && echo true || echo false) && \
		alignment::rmduplicates \
			-S ${NOrmd:=false} \
			-s ${SKIPrmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p $TMPDIR \
			-o $OUTDIR/mapped
		alignment::add4stats -r mapper && \
		alignment::bamstats \
			-S ${NOstats:=false} \
			-s ${SKIPstats:=false} \
			-r mapper \
			-t $THREADS \
			-o $OUTDIR/stats
	} || return 1

	${NOsplitreads:=true} || {
		{ 	pipeline::_slice $($sliced || ${SKIPnsplit:=false} || ${NOnsplit:=false} && echo true || echo false) && \
			alignment::splitncigar \
				-S ${NOnsplit:=false} \
				-s ${SKIPnsplit:=false} \
				-t $THREADS \
				-m $MEMORY \
				-g $GENOME \
				-r mapper \
				-c slicesinfo \
				-p $TMPDIR \
				-o $OUTDIR/mapped
		} || return 1
	}

	{	pipeline::_slice $($sliced || ${SKIPreo:=false} || ${NOreo:=false} && echo true || echo false) && \
		alignment::reorder \
			-S ${NOreo:=false} \
			-s ${SKIPreo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		pipeline::_slice $($sliced || ${SKIPlaln:=false} || ${NOlaln:=false} && echo true || echo false) && \
		alignment::leftalign \
			-S ${NOlaln:=false} \
			-s ${SKIPlaln:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped
	} || return 1

	if [[ $DBSNP ]]; then
		callvariants::vcfzip -t $THREADS -i $DBSNP || return 1
	fi

	{	pipeline::_slice $($sliced || ${SKIPbqsr:=false} || ${NObqsr:=false} && echo true || echo false) && \
		alignment::bqsr \
			-S ${NObqsr:=false} \
			-s ${SKIPbqsr:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-d "$(${NOdbsnp:-false} || echo $DBSNP)" \
			-r mapper \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		alignment::postprocess \
			-S ${NOidx:=false} \
			-s ${SKIPidx:=false} \
			-j index \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		pipeline::_slice $($sliced || ${SKIPhc:=false} || ${NOhc:=false} && echo true || echo false) && \
		callvariants::mutect \
			-S ${NOmu:=false} \
			-s ${SKIPmu:=false} \
			-t $THREADS \
			-m $MEMORY \
			-g $GENOME \
			-r mapper \
			-1 NIDX \
			-2 TIDX \
			-c slicesinfo \
			-p $TMPDIR \
			-o $OUTDIR/variants
	} || return 1

	return 0
}
