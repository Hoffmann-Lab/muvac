#! /usr/bin/env perl
use strict;
use warnings;
use v5.10;
use List::MoreUtils qw(indexes);

my $switch=1;
my $caller;
while(<>){
    my $l = $_;
    chomp $l;
    my @l = split /\s+/,$l;
    if ($l[0]=~/^CALLER=(.+)$/){
        $caller=$1;
        next;
    }
    if ($l[0]=~/^#/){
        if ($l[0]=~/^(#+)INFO=<ID=/) {
            if ($switch){
                say $1.'INFO=<ID=CALLER,Number=1,Type=String,Description="Used variant caller">';
                $switch = 0;
            } else {
                next;
            }
        } else {
            say $l;
        }
        next;
    }

    my @t = split /:/,$l[-3];
    my @vt = split /:/,$l[-1];
    my @vn = split /:/,$l[-2];
    my @tr;
    my @vtr;
    my @vnr;

    for my $k (('GT','COV','MAF','GQ','DP','AD','DP4','RSF','ASF')) {
        for (indexes { $_ eq $k } @t){
            push @tr , $t[$_];
            push @vtr , $vt[$_];
            push @vnr , $vn[$_];
        }
    }

    $l[-1] = join ':' , @vtr;
    $l[-2] = join ':' , @vnr;
    $l[-3] = join ':' , @tr;
    $l[-4] = "CALLER=$caller";

    $l = join "\t" , @l;
    say $l;
}
