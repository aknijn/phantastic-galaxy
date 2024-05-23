#!/usr/bin/env perl
#-------------------------------------
#   vcf2mst.pl 
#-------------------------------------
my $presentation=q{vcf2mst.pl         
Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree.

 - [Link to the Paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0).

When using vcf2mst please use the following citation:

SARS-CoV-2 surveillance in Italy through phylogenomic inferences based on Hamming distances derived from functional annotations of SNPs, MNPs and InDels
Adriano Di Pasquale, Nicolas Radomski, Iolanda Mangone, Paolo Calistri, Alessio Lorusso, Cesare Camma
BMC Genomics 22, 782 (2021). https://doi.org/10.1186/s12864-021-08112-0

## Synopsis

```sh
# Basic usage
vcf2mst.pl input_tsv_file out_file type_of_input [options]

# type_of_input=[vcf,gisaid,algn2pheno,nextclade,tsv,code]

# Examples

# type of inputs
vcf2mst.pl list_of_vcfiles           mst.nwk vcf
vcf2mst.pl gisaid_metadata.tsv       mst.nwk gisaid
vcf2mst.pl algn2pheno_metadata.tsv   mst.nwk algn2pheno
vcf2mst.pl nextclade_metadata.tsv    mst.nwk nextclade
vcf2mst.pl samples_vcfcodes.tsv      mst.nwk code

# profile file only. return a matrix compatible with grapetree input
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -out profile 

# filter positions
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -minmax-include 10-10000
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -minmax-exclude 10-10000

# filter positions from file
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -minmax-exclude-file intervals-file-name 10-10000

```

For additional info visit 
https://github.com/genpat-it/vcf2mst

};
#-----------------------------------
# MAIN 
#-----------------------------------
my ($f, $out, $type, $options)=@ARGV;
my %opt;

main();

sub main{
    init();
    #
    # input -> variant codes
    #
    $f=input2variantcodes($f);
    #
    # vcf2profile
    #
    $f= vcf2profile($f);


    if($opt{out} eq 'profile'){
        qx{cp $f $out};
        out( "DONE! profile file in $out \n");
    }else{
        #-----------------------------------
        # grapetree
        #-----------------------------------
        ###########profile2ids();
        _run("$grapetreeCommand $f > $out");
        ###########ids2profile();    
        out( "DONE! newick file in $out \n");
    }
}#-----------------------------------

sub input2variantcodes{ my ($f) =@_;    
    if(     $type eq 'vcf'){
                                        $f=_vcfListSnippy2Codes($f);
    }elsif( $type eq 'tsv'){
                                        $f=_tsv2Codes($f);
    }elsif( $type =~ /(nextclade|algn2pheno|gisaid)/ ){
                                        $f=_gisaidMetadata2Codes($f);
    };
    return $f;    
}#-----------------------------------

sub _vcfListSnippy2Codes{ my ($file) =@_;    
    #-----------------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from vcf files in snippy format
    # 
    # From
    #   NC_045512	241	snp	C	T	T:2201 C:0	5'UTR	+			intergenic_region n.241C>T 			
    #   NC_045512	344	snp	C	T	T:2269 C:0	mat_peptide	+	79/21290	27/7095	missense_variant c.79C>T p.Leu27Phe	GU280_gp01	ORF1ab	leader protein
    # to
    #   NC_045512      snp:C241T    241
    #   NC_045512      snp:C344T    344
    #-----------------------------------------------
    my $ar=__fileList2Array($file);
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
    my $err= '';
    foreach my $f (@{$ar}){
        if(! -e $f ){
            $err .= "# missing file: $f\n";
            next;
        }
        $cmp=$f;
        $cmp=~ s/^.*\///; # change the name of sample removing directories
        open(F, $f);
        while(<F>){
            _chomp();
            $_ =~ s/\r//g; 
            if( $_ =~ /^#/){next;}
            my $r=$_;
            # split. fields from snippy            
            # CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT
            my @ar=split(/\t/,$r); 
            #my $type=@ar[2]; 
            my ($chrom,$pos,$type,$ref,$alt)=@ar; 
            my $code='';
            if( $type eq 'snp'){
                $code= "$type:$ref$pos$alt";
            }elsif( $type eq 'mnp'){
                $code= "$type:$ref$pos$alt";
            }elsif( $type =~ 'del' ){
                $pos=$pos + 1;
                $pos2=$pos + length($ref) -2;
                $code= "$type:$pos-$pos2";
            }else{
                $code= "$type:$pos-$ref-$alt";
            }
            $out .= "$cmp\t$code\t$pos\n";
        }
    close(F);
    }
    #qx{         _vcfListSnippy2Codes.pl $f > samples_vcfcodes.csv;     };    
    my $out_file=$opttmpfile{'tmpfile-samples-vcfcodes'}; 
    open(F, ">$out_file" );
    print F $out;
    close(F);    
    return $out_file;
}#-----------------------------------


sub _gisaidMetadata2Codes{ my ($file) =@_; 
    #-------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from gisaid metadata file format
    # 
    # From
    #   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
    # to
    #   hCoV-19/Tunisia/MHT_2/2020      N_S202N     N_S202
    #   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S    NS8_L84
    # 
    #-------------------------------------
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
    open(F, $file);
    while(<F>){
        _chomp();        
        ($cmp,$vcfstring)=   __row2cmp_vcfstring($_);
        # if($_=~/^(\S+)[^\(]+\((\S+)\)/){
        if($cmp){
            #print "$cmp--- $vcfstring\n";
            @vcfs=split(/[,;]/,$vcfstring);
            my $pos;
            foreach $v (@vcfs) {
                #
                # pos calculation. compatible codes examples
                # GISAID: 
                #   NSP3_I1413L, N_S33del, Spike_ins214EPE  ->  NSP3_1413, N_33, Spike_214
                # INSA algn2pheno script:
                #   N:I1413L,    N:S33del, S:ins214EPE      ->  N:1413,    N:33, S:214
                # NEXTCLADE: 
                #   I1413L,      S33del,   ins214EPE        ->  1413,      33,   214
                #   
                $pos='-1';
                if( $v=~/^\s*([^_:]*[_:]?)(\w*\d+)/ ) {
                    $pos="$1$2";
                }
                if( $pos > 0 ){
                    $out .=  "$cmp\t$v\t$pos\n";
                }else{
                    _debug( "WARNING: row excluded in _gisaidMetadata2Codes '$_' (Not a valid position \$pos=$pos) \n" );
                }
                
            }
        }
    }
    close(F);
    
    my $out_file=$opttmpfile{'tmpfile-samples-vcfcodes'};
    open(F, ">$out_file" );
    print F $out;
    close(F);
    
    return $out_file;
}#-----------------------------------

sub _tsv2Codes{ my ($file) =@_; 
    #-------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from generic tsv format. 
    # based on options tsv-XXX 
    # 
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
    open(F, $file);
    while(<F>){
        _chomp();
        ($cmp,$vcfstring)=   __row2cmp_vcfstring($_);
        if($cmp){
            my $sep    =$opt{'tsv-mutation-sep'};
            my $regexp =$opt{'tsv-mutation-pos-regexp'};
            my $replace=$opt{'tsv-mutation-pos-replace'};
            @vcfs=split(/$sep/,$vcfstring);
            my $pos;
            foreach $v (@vcfs) {
                $pos='-1';
                $v =~ s/\r//g; # Veronica update
                if( $v=~/$regexp/ ) {                    
                    $pos=$v;
                    _debug( "DEBUG POSITION: \$pos=$pos  --> $regexp -->", 5 );
                    eval("\$pos=~ s/$regexp/$replace/");
                    _debug( "\$pos=$pos \n", 2 );
                }
                if( $pos > 0 ){
                    $out .=  "$cmp\t$v\t$pos\n";
                }else{
                    _debug( "WARNING: row excluded in _tsv2Codes '$_' (Not a valid position \$pos=$pos) \n" );
                }
            }
        }
    }
    close(F);
    #
    my $out_file=$opttmpfile{'tmpfile-samples-vcfcodes'};
    open(F, ">$out_file" );
    print F $out;
    close(F);
    return $out_file;
}#-----------------------------------



sub __row2cmp_vcfstring{  my ($s) =@_; 
    #
    #   -tsv-sample-pos pos: the position in the tsv file of column containing the sample_name (first position is 0).(default=0)
    #   -tsv-mutationslist-find string=(pos|regexp)*: the way to find the list_of_mutation_codes string in tsv file. if string=pos, -tsv-mutationslist-pos must be set.  if string=regexp, -tsv-mutationslist-regexp must be set. (default=regexp)
    #   -tsv-mutationslist-pos pos*: the position in the tsv file of column containing the list_of_mutation_codes (first position is 0) 
    #   -tsv-mutationslist-regexp string*: the regular expression used to extract the list_of_mutation_codes string. default="`\((.*)\)`"
    #   -tsv-mutation-sep char: the character used as separator between mutations on list_of_mutation_codes string. default=','
    #   -tsv-mutation-pos-regexp string*:  the regular expression used to extract the position of the mutation. default="`^(.*?[_:]?)\w*(\d+)`"
    #   -tsv-mutation-pos-replace string*:
    #
    if( $_=~ /(SAMPLECODE|TYPE|aaSubstitutions|substitutions|All\smutations)/ ){return ('','');} 
    my($cmp,$vcfstring);
    if(                                 $type eq 'gisaid'){
        if($s=~/^(\S+)[^\(]+\((\S+)\)/){
            $cmp=$1;
            $vcfstring=$2;
        }
    }elsif(                             $type eq 'algn2pheno'){
        my @a=split(/\t/,$s);
        $cmp=$a[0];
        $vcfstring=$a[9];

    }elsif(                             $type eq 'nextclade'){
        my @a=split(/\t/,$s);
        $cmp=$a[0];
        $vcfstring=$a[26];
    }elsif(                             $type eq 'tsv' ){
        if( $opt{'tsv-mutationslist-find'} eq 'pos' ){
            my $sep=$opt{'tsv-separator'};
            my $spos=$opt{'tsv-sample-pos'};
            my $mpos=$opt{'tsv-mutationslist-pos'};

            my @a=split(/$sep/,$s); 
            $cmp=$a[$spos];
            $vcfstring=$a[$mpos];
        }else{
            my $regexp=$opt{'tsv-mutationslist-regexp'};
            if($s=~/$regexp/){
                $cmp=$1;
                $vcfstring=$2;
            }
        }

    }
    return ($cmp,$vcfstring);
}#-----------------------------------

sub __fileList2Array{ my ($f) =@_;
    #-----------------------------
    # __fileList2Array
    #-----------------------------
    my @ar;
    if( $f eq '' ){
        while(<>){
            _chomp();
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
}#-----------------------------------

#############################################
#  vcf -> profile
#############################################


sub vcf2profile { my ($f)=@_;
    #-------------------------------------
    #   Hamming Distance matrix from Samples vcf codes
    #   usa:ge: 
    #   vcf2m.plst.pl samples_vcfcodes.csv > hamming_distance_matrix.tsv
    #
    #   pri ant a matrix compatible with grapetree input
    #-------------------------------------
    my $out_file=$opttmpfile{'tmpfile-samples-profiles'};

    my $c={};
    my $codes={};
    my @a; 

    open(F, $f);
    while(<F>){
        _chomp();
        $_ =~ s/\r//g; 
        # check headers
        if( $_=~ /(SAMPLECODE|TYPE)/ ){next;}
        #
        # extract code, mutations and pos
        #
        my ($SAMPLECODE,$VCFCODE,$POS) =   split(/\t/,$_);
        ###########àprint "($SAMPLECODE,$VCFCODE,$POS)\n";# =   split(/\t/,$_);
        
        if( $VCFCODE eq ''){
            ($SAMPLECODE,$VCFCODE,$POS) =   split(/,/,$_);
        }
        if( ! exists($c->{$SAMPLECODE})){
            $c->{$SAMPLECODE}={};
        }
        #
        # filter mutation  based on position 
        #
        if( _avoid_mutation_by_pos($POS) ){
            _debug("excluded: $SAMPLECODE,$VCFCODE,$POS\n",2);
            next;
        }
        # 
        # "allele" code based on variant position
        # 
        if( ! exists($codes->{$POS}) ){
            $max->{$POS}=0;
            $codes->{$POS}={};
        }
        if( ! exists($codes->{$POS}->{$VCFCODE}) ){
            my $allele_code = $max->{$POS} + 1;

            $codes->{$POS}->{$VCFCODE}=$allele_code;
            $max->{$POS}              =$allele_code;
        }
        $c->{$SAMPLECODE}->{$POS} = $codes->{$POS}->{$VCFCODE};
    }
    close(F);
    #---------------------------
    # out
    #---------------------------
    my $out="";
    #
    # Header
    #
    $out .= "#FILE";
    foreach $cod (sort(keys(%$codes))){    $out .= "\t$cod";}
    $out .= "\n";
    #
    # "allele" codes
    #
    foreach $cmp (sort(keys(%$c))){
        $out .= "$cmp";
        foreach $cod (sort(keys(%$codes))){
            my $allele=$c->{$cmp}->{$cod};
            my $val=( $allele )?
                                $allele + 1 :
                                1;
            $out .= "\t$val";
        }
        $out .= "\n";
    }
    #
    # print
    #
    open(F, ">$out_file" );
    print F $out;
    close(F);    
    return $out_file;
}#-----------------------------------

sub _avoid_mutation_by_pos{ my ($pos)=@_;
    #
    # TODO: manage -minmax-include and -minmax-exclude
    #    
    if(_containsElements($hs_minmax_exclude)){
        if(__check_minmax($pos,$hs_minmax_exclude)){
            _debug( "$pos excluded\n",3);
            return 1;                                   #1 avoid=exclude
        }
    }
    if(_containsElements($hs_minmax_include)){
        if(__check_minmax($pos,$hs_minmax_include)){
            _debug( "$pos included\n",1);
            return 0;                                   #0 not avoid=include
        }else{
            _debug( "$pos excluded\n",1);
            return 1;                                   #avoid=exclude
        }
    }
    return 0; #not avoid=include
}#-----------------------------------

sub __check_minmax{ my ($pos, $hs)=@_;
    #use Data::Dumper; print Dumper($hs);  
    foreach my $h (values(%{$hs})){
        if( 
            $pos >= $h->{min} && 
            $pos <= $h->{max}
        ){
            _debug( "DEBUG(l3) pos=$pos in interval[$h->{min},$h->{max}]. ",5);
            return 1 # included in the interval
        }
    }
    return '';
}#-----------------------------------



sub init{
    #--------------------------------------
    #   INIT
    #--------------------------------------
    #
    # input filename
    # 
    if( ! $f ){
        print $presentation;
        exit;          
    }    
    #
    # type_of_input
    # 
    if( $type =~ /^(vcf|gisaid|algn2pheno|nextclade|tsv|code)$/ ){
        _debug( "type_of_input: $type \n");
    }else{
        print "ERROR: MODE NOT RECOGNIZED: $type\n";
        exit;
    }
    #
    # output filename
    # 
    if( ! $out ) {
        $out='/tmp/mst.nwk';
    }
    #
    # options
    # 
    _init_set_option_defaults();
    my $initial_options=  join(' ',@ARGV);
    $initial_options   =~ s/^.*?-/-/;    
    $options=$initial_options;
    while( $options =~ /^\s*-(\S+)\s+(\S+)(.*)$/ ){
        #print "---options=$options\n";
        $opt{$1}=$2;
        $options=$3;
    }
    if( $initial_options ne $options){
        _debug( "find options:\n");
        foreach my $k (keys(%opt)){
            _debug( "     -$k=$opt{$k}\n");
            if(   $k eq 'out' ){
                if($opt{$k} !~ /profile/){
                    _debug( qq{ERR: value "$opt{$k}" not known for option "$k"
                    example: -out profile\n});
                    exit;
                }
            }elsif($k =~ /^(minmax-include|minmax-exclude)$/ ){
                if($opt{$k} !~ /^([^-]+-[^-]+,?)+$/){
                    _debug( qq{ERR: value "$opt{$k}" not correct for option "$k"
                    example: -minmax-include 0-1000
                    example: -minmax-exclude 0-1000,1200-12111
                    \n});
                    exit;
                }
                # check for wrong use minmax while the idea was to use file-minmax
                if(-e $opt{$k} ){
                    _debug( qq{ERR: value "$opt{$k}" is a filename. not correct for option "$k"
                    Did you intend to use -file-$k ?
                    \n});
                    exit;
                }
            }elsif($k =~ /^(file-minmax-include|file-minmax-exclude)$/ ){
                my $file=$opt{$k};
                if(! -e $opt{$k} ){
                    _debug( qq{ERR: "$opt{$k}" is not a file. option "$k"
                    example: -file-minmax-include           file_of_intervals
                    example: -file-minmax-exclude   file_of_intervals
                    \n});
                    exit;
                }
            }elsif($k =~ /^(tsv-separator|tsv-sample-pos|tsv-mutationslist-find|tsv-mutationslist-pos|tsv-mutationslist-regexp|tsv-mutation-sep|tsv-mutation-pos-regexp|tsv-mutation-pos-replace)$/){
                #ok tsv options
            }elsif($k =~ /^debug$/){
                if($opt{$k} !~ /^\d+$/ ){
                    print qq{ERR: An integer (debug level) must specified for option "$k". Value "$opt{$k}" is not correct 
                    example: -debug 1
                    \n};
                    exit;
                }
            }elsif($k =~ /^grapetree-bin$/){
                #ok                
            }
        }
    }
    if($options && ($options ne $initial_options)){
        print "options not recognized\n$options\n";
        exit;
    }
    #
    #
    #
    _init_include_exclude_pos();    
    #
    # grapetree-bin
    # 
    $grapetreeCommand=_init_grapetree();
}#-----------------------------------


sub _init_include_exclude_pos {
    #
    # 
    #
    $hs_minmax_include=__read_include_exclude_pos('file-minmax-include','minmax-include');
    $hs_minmax_exclude=__read_include_exclude_pos('file-minmax-exclude','minmax-exclude');
    
    _debug( 'DEBUG list minmax_include='.join(',',keys(%{$hs_minmax_include})) . "\n", 1);
    _debug( 'DEBUG list minmax_exclude='.join(',',keys(%{$hs_minmax_exclude})) . "\n", 1);
    
}#-----------------------------------

sub __read_include_exclude_pos{ my ($f,$s)=@_; 
    # f= file-minmax-include or file-minmax-exclude
    # s= minmax-include      or minmax-exclude
    # 1. read from file (file-minmax-xx). 
    # 2. add to string (minmax-xx)
    # 3. create hash of intervals
    if($opt{$f}){
        open(F,$opt{$f});
        my $strFromFile='';
        while(<F>){
            _chomp();
            $_ =~ s/\r//g; 
            # empty or a line commented (#)
            if( $_ =~/^\s*#/ || $_ =~ /^\s*$/ ) { next; }
            #remove spaces
            $_ =~ s/[\s\t\r\n]+//g;
            #add with comma
            $strFromFile.= "$_,";
        }
        close(F);
        if( $strFromFile  ){
            $opt{$s} = ( $opt{$s} )? 
                        "$opt{$s},$strFromFile" : 
                        $strFromFile;
        }
    }
    if($opt{$s}){
        return __get_minmax($opt{$s});
    }    
    return {};
}#-----------------------------------

sub __get_minmax{ my ($minmax_list_string)=@_;
    my @ar_minmax_list_string=split(/,/,$minmax_list_string);
    my $hs_out={};
    foreach my $mm (@ar_minmax_list_string){
        if($mm =~ /^\s*$/){next;}
        my ($min,$max)=split(/-/,$mm);
        $min =~ s/[\s\t\r\n]*//g;
        $max =~ s/[\s\t\r\n]*//g;
        $hs_out->{$mm}={
            min => $min, max => $max 
        };
    }
    #use Data::Dumper;print Dumper($hs_out);    
    return $hs_out;
}#-----------------------------------

sub _init_grapetree {
    #
    # set grapetree bin 
    # based on -out options and GRAPETREE_EXEC env variable
    # example of $ENV{GRAPETREE_EXEC}='docker run  --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p';
    # 
    if($opt{out} eq 'profile'){
        _debug( "_init_grapetree: no need of grapetree for options -out profile \n");
        return '';
        # no need of grapetree for options -out profile
    }
    if( ! $ENV{GRAPETREE_EXEC} ){
        _debug( "_init_grapetree: env GRAPETREE_EXEC is NOT set. grapetree will be used \n");
        return $opt{'grapetree-bin'};
    }else{
        _debug( "_init_grapetree: env GRAPETREE_EXEC is set: '$ENV{GRAPETREE_EXEC}'\n");
        return $ENV{GRAPETREE_EXEC};
    }
}#-----------------------------------

sub _init_set_option_defaults {
    #
    # tsv options
    #
    $opt{'tsv-separator'}="\t";
    $opt{'tsv-sample-pos'}=0;
    $opt{'tsv-mutationslist-find'}='regexp';
    $opt{'tsv-mutationslist-pos'}=0;
    $opt{'tsv-mutationslist-regexp'}='\((.*)\)';
    $opt{'tsv-mutation-sep'}=',';
    $opt{'tsv-mutation-pos-regexp'} ='^(.*?[_:]?)\w*(\d+)';
    $opt{'tsv-mutation-pos-replace'}='$1$2';
    #
    $opt{'debug'}=0;
    #
    # grapetree bin
    #
    $opt{'grapetree-bin'}='grapetree -p ';
    #
    # tmpfiles
    #
    my $t=time;
    my $name='tmpfile-samples-vcfcodes';
    $opttmpfile    {$name}="/tmp/samples-vcfcodes-$t.tsv";
    $opttmpfileDesc{$name}='full list of recognized mutation codes (not filtered by position)';
    #
    my $name='tmpfile-samples-profiles';
    $opttmpfile    {$name}="/tmp/samples-profiles-$t.tsv"; #ex /tmp/hamming_distance_matrix
    $opttmpfileDesc{$name}='mutations codes in profile format (possibly filtered by position)'
}#-----------------------------------



#-----------------------------------
# OUT
#-----------------------------------
sub out{ my ($s) =@_;    
    print "tmp files in  \n";
    foreach my $k (keys(%opttmpfile)){
        print "  $opttmpfile{$k} : $opttmpfileDesc{$k}\n";
    }
    print $s;
}#-----------------------------------


#-----------------------------------
# BASIC UTILS
#-----------------------------------
sub _run{ my ($s) =@_;    
    #
    #   run a bash command taking start stop time
    #
    print "$s";
    print "\nstart: " . qx{date};
    qx{$s};
    print "\nstop: " . qx{date};
}#-----------------------------------

sub _debug{ my ($s, $level) =@_;    
    #
    #   print string if debug is on (option -debug on)
    #
    if( !$level ){$level=0;}
    my $debug_level=$opt{debug};
    if( $level < $opt{debug}){
        print "$s";
    }
}#-----------------------------------

sub _chomp{
    #
    #   remove \r and \n from $_ (base chomp statement remove just \n)
    #
    chomp;
    $_ =~ s/\r//g; 
}#-----------------------------------

sub _containsElements{my ($s) =@_;    
    #
    #   return true if the reference to hash is not empty
    #
    return %{$s};
}#-----------------------------------



<<_________COMMENT_____________;

/tmp/samples_vcfcod
/tmpes
/tmp/hamming_distance_matr
/tmpix
/tmp/m.nwkst.nwk
_________COMMENT_____________




