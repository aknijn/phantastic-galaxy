#!/usr/bin/env perl
#
#-------------------------------------
#   vcf2mst.pl 
#  
#   Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree
#   usage 1: 
#   vcf2mst.pl samples_vcfcodes.csv mst.nwk 
#
#   usage 2: 
#   vcf2mst.pl list_of_vcfiles mst.nwk  vcf
#
#-------------------------------------
#
my ($f, $out, $type)=@ARGV;

if($type eq 'vcf'){
    qx{
        vcflist2codes.pl $f > samples_vcfcodes.csv;
    };
    $f='samples_vcfcodes.csv';
};

# esecuzione 
# 
run("perl vcf2ham.pl $f > /tmp/hdmatrix.tsv");
run("docker run --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p /tmp/hdmatrix.tsv > /tmp/mst.nwk ");

if( $out ) {
    qx{mv /tmp/mst.nwk  $out};
}else{
    $out='/tmp/mst.nwk';
};

print "newick file in $out \n";

#-----------------------------------
# FUNZIONI
#-----------------------------------

sub run{ my ($s) =@_;    
    print "$s";
    print "\nstart: " . qx{date};
    qx{$s};
    print "\nstop: " . qx{date};
}

#sudo docker run --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p /tmp/ham.tsv > /tmp/a.nwk
#  graptree -p hdmatrix.tsv > mst.nwk 

# qx{    perl vcf2ham.pl $f > /tmp/hdmatrix.tsv;
#    sudo docker run --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p /tmp/hdmatrix.tsv > /tmp/mst.nwk  };

