#!/usr/bin/env perl

my ($f)=@ARGV;

if(! $f){die(q{vcflist2codes.pl 

    From a list of VCF files (snippy output format) to TSV of SAMPLECODE,VCFCODE
    
    usage 1: 
    vcflist2codes.pl list_of_vcf_files  > samples_vcfcodes.tsv

    list_of_vcf_files can be 

    1. a filename F:
        F contains a list of filenames divided by newline (e.g., ls dir > list_of_vcf_files). 

    2. a directory D:
        A `find D` will be executed for input filenames


    Each file in the list must have the vcf format like this (e.g., from snippy output):
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  DS11243446_vdsnippy_NC045512
    NC_045  241     .       C       T       42316.5 .       AB=0;AO=1409;DP=1411;QA=47794;QR=36;RO=1;TYPE=snp;ANN=T|intergenic_region|MODIFIER|CHR_START-GU280_gp01|CHR_START-GENE_GU280_gp01|intergenic_region|CHR_START-GENE_GU280_gp01|||n.241C>T||||||  GT:DP:RO:QR:AO:QA:GL    1/1:1411:1:36:1409:47794:-47794

    The filename without directory name will be used for SAMPLECODE

});
}

#=============================
# MAIN
#=============================
my $ar=array_of_vcf_files($f);


foreach my $f (@{$ar}){
    if(! -e $f ){
        print 
    }

    $cmp=$f;
    $cmp=~ s/^.*\///; # change the name of sample removing directories
    print "SAMPLECODE\tVCFCODE\n";

    open(F, $f);
    while(<F>){
 	chomp;
     	    if( $_ =~ /^#/){next;}

        my $r=$_;
        my @ar=split(/\t/,$r);
        my $type=''; my $code='';
        if( $r =~ /;TYPE=(\w+);/){
            $type=$1;
        }

        if( $type eq 'snp'){
            $code= "$type:@ar[3]@ar[1]@ar[4]";
        }elsif( $type eq 'mnp'){
            $code= "$type:@ar[3]@ar[1]@ar[4]";
        }elsif( $type =~ 'del' ){
            $pos=@ar[1] + 1;
            $pos2=$pos + length(@ar[3]) -2;
            $code= "$type:$pos-$pos2";
        }else{
            $code= "$type:@ar[1]-@ar[3]-@ar[4]";
        }
        print "$cmp\t$code\n";
    }
    close(F);
}

#-----------------------------
# array of files
#-----------------------------
sub array_of_vcf_files{ my ($f) =@_;
    my @ar;
    if( $f eq '' ){
        while(<>){
            chomp;
            push(@ar, $_);
        }
    }elsif( -e $f ){
        my $ris='';
        if( -d $f ){
            $ris=qx{find $f};
        }else{
            $ris=qx{cat $f};
        }
        @ar=split(/\n/,$ris);
	# use Data::Dumper;         print Dumper(\@ar);
    }else{
        die("file/dir not present: $f ")
    }
    return \@ar;
};

