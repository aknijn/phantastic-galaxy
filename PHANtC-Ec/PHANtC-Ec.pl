#!/usr/bin/env perl
## A wrapper script to call chewBBACA and MentaLiST (only tree part)
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
    $serotype_file,
    $phantcec_allele,
    $phantcec_json,
    $phantcec_tree,
    $phantcec_dm) = @ARGV;

# Run program
my $chewiedir = '/gfs/data-flow/Chewie-NS';
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $sampleGenesMapped = "0";
my $cgLociNumber = 0;
my $permille_loci = 0;
prepareEnvironment($input1,$input1_name,"input_dir");
runChewBBACA();
collectStatistics();
collectOutput();
runMentaLiST();
exit(0);

# Run chewBBACA
sub runChewBBACA {
    my $allelecalldir = "$scriptdir/allelecall";
    my $utilsdir = "$scriptdir/utils";
    my $newpath = "PATH=$ENV{PATH}:$allelecalldir:$utilsdir";
    my $python = "chewBBACA.py AlleleCall -o output_dir -i input_dir --cpu 4 --hash-profiles crc32 --no-inferred --bsr 0.6 --ptf $chewiedir/prodigal_training_files/Escherichia_coli.trn -g $chewiedir/ecoli/ecoli_INNUENDO_wgMLST_ORIG/ --gl $chewiedir/ecoli/Ecoli_cgMLST_ns_ids.txt";
    my $result = system("$newpath; $python");
    return 0;
}

# Run runMentaLiST, only the tree part
sub runMentaLiST {
    createAllelesFile();
    open my $of, '>', $phantcec_json or die "Cannot open json: $!";
    print $of "{\"core_genome_schema_size\": $cgLociNumber, \"sample_genes_mapped\": $sampleGenesMapped}";
    close $of;
    # calc distance matrix from database
    my $result = system("python $scriptdir/scripts/mlst_hash_stretch_distance.py -i cgMLST.tmp -o $phantcec_dm");
    # calc tree from distance matrix
    system("python $scriptdir/scripts/mentalist_tree $phantcec_dm > $phantcec_tree");
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
    # remove INF- from newly inferred alleles, substitute - with 0 and remove .fasta from the filename
    $allele_line =~ s/INF-//ig;
    $allele_line =~ s/-/0/ig;
    $allele_line =~ s/.fasta//ig;
    close $if_in;

    my $sql_insert = "insert into mlst_ecoli (sample_code) values (?)";
    my $sql_update = "update mlst_ecoli set permille_loci=?, allele_strain=? where sample_code=?";
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
    return 0;
}

# Collect output with statistics, save the number of genes mapped, the total number of loci and the relative number of mapped genes
sub collectStatistics{
    move("output_dir/results_statistics.tsv", $phantcec_allele) ;
    open(my $if_st, '<', $phantcec_allele) or die "Could not read from results_statistics.tsv, program halting.";
    <$if_st>;
    my $statistics_line = <$if_st>;
    chomp $statistics_line;
    my @keys = split( /\t/, $statistics_line );
    $sampleGenesMapped = int($keys[1]) + int($keys[2]);
    $cgLociNumber = int($keys[1]) + int($keys[2]) + int($keys[3]) + int($keys[4]) + int($keys[5]) + int($keys[6]) + int($keys[7]) + int($keys[8]) + int($keys[9]) + int($keys[10]) + int($keys[11]);
    $permille_loci = int((int($sampleGenesMapped)/$cgLociNumber)*1000 + 0.5);
    close $if_st;      
    return 0;
}

# Obtain allele profiles from db with allele_coverage>79,9% and write to temp file
sub createAllelesFile {
    open(my $if, '<', $serotype_file) or die "Could not read from serotype_file, program halting.";
    my $serotype = <$if>;
    chomp $serotype;
    close $if;
    my $sql;
    # If the serogroup is known, retrieve allele profiles only from those samples, if O? then from all
    if ($serotype eq 'O?') {
        $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
               select allele_strain from mlst_ecoli where permille_loci>799";
    }
    else
    {
        $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
               select allele_strain from mlst_ecoli where sample_code='$input1_name' union
               select allele_strain from v_mlst_ecoli where permille_loci>799 and Antigen_O = '$serotype'";
    }
    # connect to MySQL database
    my %attr2 = ( PrintError=>0, RaiseError=>1);
    my $dbh2 = DBI->connect($dsn,$user,$pwd,\%attr2);
    my $sth2 = $dbh2->prepare($sql);
    $sth2->execute();

    open my $of, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
    while(my @row = $sth2->fetchrow_array()){
      print $of "$row[0]\n"
    }       
    close $of;

    $sth2->finish();
    # disconnect from the MySQL database
    $dbh2->disconnect();
    return 0;
}