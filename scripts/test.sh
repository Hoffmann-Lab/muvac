#! /usr/bin/env bash

THREADS=$(grep -cF processor /proc/cpuinfo)
[[ $MUVAC ]] || {
	echo -n "where is muvac? [/path/to/install/dir]: "
	read -r p
	export MUVAC=$(readlink -e ${p/#*~/$HOME})
}

### preprocessing
prep(){
	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-o /ssd/tmp/muvac_test_se

	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-2 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
		-o /ssd/tmp/muvac_test_pe
}
echo -n "test preprocessing? [y|n]: "
read -r yn
[[ $yn == "y" ]] && prep

### mapping
map(){
	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-g /misc/paras/data/genomes/GRCh37.p13/GRCh37.p13.fa \
		-1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-2 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
		-o /ssd/tmp/muvac_test_pe \
		-tmp /ssd/tmp/muvac_test_pe \
		-resume hc \
		-skip md5
	exit

	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-g /misc/paras/data/genomes/GRCh37.p13/GRCh37.p13.fa \
		-1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-o /ssd/tmp/muvac_test_se \
		-tmp /ssd/tmp/muvac_test_se \
		-resume hc \
		-skip md5,slice
	exit



}
echo -n "test mapping? [y|n]: "
read -r yn
[[ $yn == "y" ]] && map

