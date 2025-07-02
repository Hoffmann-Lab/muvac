#! /usr/bin/env bash
# (c) Konstantin Riege

function compile::bashbone(){
	local insdir threads cfg version bashboneversion src="$(dirname "$(dirname "$(readlink -e "$0")")")"
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$src/lib/version.sh"
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	rm -f "$insdir/bashbone-$version/scripts/"+(setup|test).sh
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	bashboneversion=$version
	src="$(dirname "$src")"

	commander::printinfo "installing muvac"
	source "$src/lib/version.sh"
	rm -rf "$insdir/muvac-$version"
	mkdir -p "$insdir/muvac-$version"
	cp -r "$src"/!(bashbone|setup*) "$insdir/muvac-$version"
	rm -f "$insdir/muvac-$version/scripts/"+(setup|test).sh
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/muvac-$version" "$insdir/latest/muvac"
	ln -sfn "$insdir/bashbone-$bashboneversion" "$insdir/muvac-$version/bashbone"

	return 0
}

function compile::muvac(){
	compile::bashbone "$@"

	return 0
}