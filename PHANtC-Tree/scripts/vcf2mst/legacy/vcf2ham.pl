#!/usr/bin/env perl
#
#-------------------------------------
#   vcf2ham.pl 
#  
#   Hamming Distance matrix from Samples vcf codes
#   usage: 
#   vcf2mst.pl samples_vcfcodes.csv > hamming_distance_matrix.tsv
#
#   print a matrix compatible with grapetree input
#-------------------------------------
my ($f)=@ARGV;

my $c={};
my $codes={};
my @a; 


open(F, $f);
while(<F>){
    chomp;
    if( $_=~ /SAMPLECODE/ ){next;}

    my ($SAMPLECODE,$VCFCODE) =   split(/\t/,$_);
    if( $VCFCODE eq ''){
    	($SAMPLECODE,$VCFCODE) =   split(/,/,$_);
    }
    
    if( ! exists($c->{$SAMPLECODE})){
            $c->{$SAMPLECODE}={};
    }

    $c->{$SAMPLECODE}->{$VCFCODE}=1;
    $codes->{$VCFCODE}++;
}

close(F);
#--------------------------------------

#--------------------------------------
print "#FILE\t";
foreach $cod (sort(keys(%$codes))){    print "$cod\t";}
print "\n";

foreach $cmp (sort(keys(%$c))){
    print "$cmp\t";
    foreach $cod (sort(keys(%$codes))){    
        my $val=( $c->{$cmp}->{$cod} )?1:2;
        print"$val\t";
    }
    print "\n";
}
#--------------------------------------
