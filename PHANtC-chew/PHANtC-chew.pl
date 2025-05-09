#!/usr/bin/env perl
## A wrapper script to call chewBBACA
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use DBI;
use Config::Simple;

# Parse arguments
my ($input1,
    $input1_name,
	$species,
    $phantcchew_alleleprofile,
    $phantcchew_allele,
    $phantcchew_json) = @ARGV;

# Run program
my $chewiedir = '/disk2/Chewie-NS';
my $python = "";
if ($species eq "Listeria monocytogenes") {
    my $ptfdir = "$chewiedir/prodigal_training_files/Listeria_monocytogenes.trn";
	my $cgmlstdir = "$chewiedir/listeria/lmonocytogenes_Pasteur_cgMLST_ORIG/";
	$python = "chewBBACA.py AlleleCall -o output_dir -i input_dir --cpu 4 --hash-profiles crc32 --no-inferred --bsr 0.6 --ptf $ptfdir -g $cgmlstdir";
}
elsif ($species eq "Escherichia coli") {
	my $ptfdir = "$chewiedir/prodigal_training_files/Escherichia_coli.trn";
	my $wgmlstdir = "$chewiedir/ecoli/ecoli_INNUENDO_wgMLST_ORIG/";
	my $cgmlstfilter = "$chewiedir/ecoli/Ecoli_cgMLST_ns_ids.txt";
    $python = "chewBBACA.py AlleleCall -o output_dir -i input_dir --cpu 4 --hash-profiles crc32 --no-inferred --bsr 0.6 --ptf $ptfdir -g $wgmlstdir --gl $cgmlstfilter";
}

my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $permille_loci = 0;
prepareEnvironment($input1,$input1_name,"input_dir");
runChewBBACA();
collectStatistics();
collectOutput();
exit(0);

# Run chewBBACA
sub runChewBBACA {
    my $allelecalldir = "$scriptdir/allelecall";
    my $utilsdir = "$scriptdir/utils";
    my $newpath = "PATH=$ENV{PATH}:$allelecalldir:$utilsdir";
    my $result = system("$newpath; $python");
    return 0;
}

# Run prepareEnvironment, create the directory $indir with symlink to the file in $inlist
sub prepareEnvironment {
    my ($inlist, $inlist_name, $indir) = @_;
    if ($inlist ne "NULL") {
      mkdir($indir);
      if ($inlist_name ne "NULL") {
        symlink($inlist, $indir . "/" . $inlist_name .  ".fasta");
      }
    }
    return 0;
}

# Collect output to database
sub collectOutput{
    open(my $if_in, '<', "output_dir/results_alleles_hashed.tsv") or die "Could not read from results_alleles_hashed.tsv, program halting.";
    <$if_in>;
    my $allele_line = <$if_in>;
    chomp $allele_line;
    # remove INF- from newly inferred alleles, substitute - with 0 and remove .fasta from the filename to recover sample_code
    $allele_line =~ s/INF-//ig;
    $allele_line =~ s/-/0/ig;
    $allele_line =~ s/.fasta//ig;
    close $if_in;
    open my $of_ap, '>', $phantcchew_alleleprofile or die "Cannot open alleleprofile: $!";
    print $of_ap $allele_line;
    close $of_ap;
	# only insert samples that are not (proficiency) test nor accreditation samples
	my $sample_type = substr $input1_name, 0, 2;
    if ($sample_type ne "PT" && $sample_type ne "TT" && $sample_type ne "AC") {
        my $sql_insert = "";
        my $sql_update = "";
        if ($species eq "Listeria monocytogenes") {
            $sql_insert = "insert into mlst_listeria (sample_code) values (?)";
            $sql_update = "update mlst_listeria set permille_loci=?, allele_strain=? where sample_code=?";
        }
        elsif ($species eq "Escherichia coli") {
            $sql_insert = "insert into mlst_ecoli (sample_code) values (?)";
            $sql_update = "update mlst_ecoli set permille_loci=?, allele_strain=? where sample_code=?";
        }
        # connect to MySQL database
        my %attr = ( PrintError=>0, RaiseError=>0);
        my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
        my $sth_insert = $dbh->prepare($sql_insert);
        $sth_insert->execute($input1_name);
        my $sth_update = $dbh->prepare($sql_update);
        $sth_update->bind_param( 1, $permille_loci );
        $sth_update->bind_param( 2, $allele_line );
        $sth_update->bind_param( 3, $input1_name );
        $sth_update->execute();
        $sth_insert->finish();
        $sth_update->finish();
        # disconnect from the MySQL database
        $dbh->disconnect();
    }
    return 0;
}

# Collect output with statistics, save the number of genes mapped, the total number of loci and the relative number of mapped genes
sub collectStatistics{
    move("output_dir/results_statistics.tsv", $phantcchew_allele) ;
    open(my $if_st, '<', $phantcchew_allele) or die "Could not read from results_statistics.tsv, program halting.";
    <$if_st>;
    my $statistics_line = <$if_st>;
    chomp $statistics_line;
    my @keys = split( /\t/, $statistics_line );
    my $sampleGenesMapped = int($keys[1]) + int($keys[2]);
    my $cgLociNumber = int($keys[1]) + int($keys[2]) + int($keys[3]) + int($keys[4]) + int($keys[5]) + int($keys[6]) + int($keys[7]) + int($keys[8]) + int($keys[9]) + int($keys[10]) + int($keys[11]);
    $permille_loci = int((int($sampleGenesMapped)/$cgLociNumber)*1000 + 0.5);
    close $if_st;      
    open my $of, '>', $phantcchew_json or die "Cannot open json: $!";
    print $of "{\"core_genome_schema_size\": $cgLociNumber, \"sample_genes_mapped\": $sampleGenesMapped}";
    close $of;
    return 0;
}
