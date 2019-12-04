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
		-o /ssd/tmp/muvac_germ_se

	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-2 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
		-o /ssd/tmp/muvac_germ_pe

	$MUVAC/latest/muvac/muvac.sh \
		-v 2 \
		-t $THREADS \
		-n1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
		-n2 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
		-o /ssd/tmp/muvac_som_se
}
echo -n "test preprocessing? [y|n]: "
read -r yn
[[ $yn == "y" ]] && prep

pon(){
    $MUVAC/latest/muvac/muvac.sh \
        -mem 20000 \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
        -o /ssd/tmp/muvac_som_se \
        -no-stats \
        -pon \
        -no-pondb \
        -resume pon \
        -skip md5
    $MUVAC/latest/muvac/muvac.sh \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -1 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
        -o /ssd/tmp/muvac_som_se \
        -pon \
        -no-pondb \
        -resume pon \
        -skip md5
}
echo -n "test pon? [y|n]: "
read -r yn
[[ $yn == "y" ]] && pon

pondb(){
    $MUVAC/latest/muvac/muvac.sh \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -1 $MUVAC/latest/muvac/data/test.R1.fastq.gz,$MUVAC/latest/muvac/data/test.R2.fastq.gz \
        -o /ssd/tmp/muvac_som_se \
        -pon \
        -resume pondb \
        -skip md5
}
echo -n "test pondb? [y|n]: "
read -r yn
[[ $yn == "y" ]] && pondb

call(){
	$MUVAC/latest/muvac/muvac.sh \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
        -o /ssd/tmp/muvac_germ_se \
        -resume hc \
        -skip md5

	$MUVAC/latest/muvac/muvac.sh \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
        -2 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
        -o /ssd/tmp/muvac_germ_pe \
        -resume hc \
        -skip md5

    $MUVAC/latest/muvac/muvac.sh \
        -v 2 \
        -t $THREADS \
        -g /misc/paras/data/genomes/GRCh38.p12/GRCh38.p12.fa \
        -n1 $MUVAC/latest/muvac/data/test.R1.fastq.gz \
        -t1 $MUVAC/latest/muvac/data/test.R2.fastq.gz \
        -o /ssd/tmp/muvac_som_se \
        -resume mu \
        -mypon \
        -skip md5
}
echo -n "test calling? [y|n]: "
read -r yn
[[ $yn == "y" ]] && call
