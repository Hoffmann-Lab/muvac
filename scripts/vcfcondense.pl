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
        if ($l[0]=~/^(#+)CHROM/) {
            say "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$caller";
        } elsif ($l[0]=~/^(#+)INFO=<ID=/) {
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

    my @t = split /:/,$l[-2];
    my @v = split /:/,$l[-1];
    my @tr;
    my @vr;

    for my $k (('GT','COV','MAF','GQ','DP','AD','DP2')) {
        for (indexes { $_ eq $k } @t){
            push @tr , $t[$_];
            push @vr , $v[$_];
        }
    }

    $l[-1] = join ':' , @vr;
    $l[-2] = join ':' , @tr;
    $l[-3] = "CALLER=$caller";

    $l = join "\t" , @l;
    say $l;
}
