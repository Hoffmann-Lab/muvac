#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::muvac -i "$insdir" -t $threads
	compile::conda -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads
	compile::java -i "$insdir" -t $threads
	compile::trimmomatic -i "$insdir" -t $threads
	compile::sortmerna -i "$insdir" -t $threads
	compile::segemehl -i "$insdir" -t $threads

	return 0
}

compile::muvac() {
	local insdir threads version bashboneversion src=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
	commander::printinfo "installing muvac"
	compile::_parse -r insdir -s threads "$@"
	#source $src/bashbone/lib/version.sh
	source $insdir/latest/bashbone/lib/version.sh
	bashboneversion=$version
	source $src/lib/version.sh
	rm -rf $insdir/muvac-$version
	mkdir -p $insdir/muvac-$version
	cp -r $src/!(bashbone|setup*) $insdir/muvac-$version
	mkdir -p $insdir/latest
	ln -sfn $insdir/muvac-$version $insdir/latest/muvac
	ln -sfn $insdir/bashbone-$bashboneversion $insdir/muvac-$version/bashbone

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::muvac -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads -u true

	return 0
}

compile::conda_tools() {
	local insdir threads upgrade=false url version tool n bin doclean=false
	declare -A envs

	compile::_parse -r insdir -s threads -c upgrade "$@"
	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r tool; do
		envs[$tool]=true
	done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

	# python 3 envs
	# ensure star compatibility with CTAT genome index and STAR-fusion
	for tool in fastqc cutadapt rcorrector star=2.7.2b bwa picard bamutil gatk4 freebayes varscan igv; do
		n=${tool//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} && continue
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		# link commonly used base binaries into env
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	done
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

	tool=vardict
	n=${tool//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool vardict-java
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=snpeff
	n=${tool//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool snpsift
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	# python 2 envs
	tool=platypus-variant
	n=platypus
	$upgrade && ${envs[$n]:=false} || {
		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=2
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	$doclean && {
		commander::printinfo "conda clean up"
		conda clean -y -a
	}
	conda deactivate

	return 0
}
