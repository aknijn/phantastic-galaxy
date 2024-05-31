#!/usr/bin/env perl
#
#-------------------------------------
# da
#   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
# a
#   hCoV-19/Tunisia/MHT_2/2020      N_S202N
#   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S
#
# usage 1: 
#   gisad2vcf.pl  samples_vcfcodes.csv 
#
#-------------------------------------
#

while(<>){
    #print $_;
    if($_=~/^(\S+)\s+\((\S+)\)/){
        $cmp=$1;
        $vcfstring=$2;
        #print "$cmp--- $vcfstring\n";
        @vcfs=split(/,/,$vcfstring);
        foreach $v (@vcfs) {
            print "$cmp\t$v\n";
        }
    }
}
