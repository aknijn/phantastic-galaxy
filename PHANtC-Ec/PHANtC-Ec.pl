#!/usr/bin/env perl
## A wrapper script to call chewBBACA and MentaLiST
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
    my $python = "python $scriptdir/PHANtC-Ec.py -o output_dir -i input_dir --cpu 4 --bsr 0.6 --ptf $scriptdir/TrainingFiles4Prodigal/trained_eColi.trn -g $scriptdir/../../INNUENDO/schema_chewBBACA_cgMLST_V4/_schema.txt";
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
        symlink($inlist, $indir . "/" . $inlist_name);
      }
    }
    return 0;
}

# Collect output to database
sub collectOutput{
    my @alleles = glob "output_dir/results_*/results_alleles.tsv";
    if (@alleles == 1) {
      open(my $if_in, $alleles[0]) or die "Could not read from results_alleles.tsv, program halting.";
      <$if_in>;
      my $allele_line = <$if_in>;
      chomp $allele_line;
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
    }
    return 0;
}

# Collect output with statistics, save the number of genes mapped, the total number of loci and the relative number of mapped genes
sub collectStatistics{
	my @statistics = glob "output_dir/results_*/results_statistics.tsv";
    if (@statistics == 1) { 
	  move($statistics[0], $phantcec_allele) ;
      open(my $if_st, '<', $phantcec_allele) or die "Could not read from results_statistics.tsv, program halting.";
      my $lastline;
      $lastline = $_, while (<$if_st>);
	  chomp $lastline;
	  my @keys = split( /\t/, $lastline );
	  $sampleGenesMapped = $keys[1];
	  $cgLociNumber = int($keys[1]) + int($keys[2]) + int($keys[3]) + int($keys[4]) + int($keys[5]) + int($keys[6]) + int($keys[7]);
	  $permille_loci = int((int($sampleGenesMapped)/$cgLociNumber)*1000 + 0.5);
      close $if_st;	  
	}
    return 0;
}

# Obtain alleles from db and write to temp file
sub createAllelesFile {
    open(my $if, $serotype_file) or die "Could not read from serotype_file, program halting.";
    my $serotype = <$if>;
    chomp $serotype;
    close $if;
	my $sql;
	if ($serotype == 'O?') {
		$sql = "select mlst_ecoli.allele_strain from mlst_ecoli where sample_code='FILE' union
               select mlst_ecoli.allele_strain from mlst_ecoli where permille_loci>799";
	}
	else
	{
		$sql = "select mlst_ecoli.allele_strain from mlst_ecoli where sample_code='FILE' union
               select mlst_ecoli.allele_strain from mlst_ecoli where sample_code='$input1_name' union
               select mlst_ecoli.allele_strain from sample_metadata_entry AS sme_s
                 inner join metadata_field AS mf_s on sme_s.metadata_KEY = mf_s.id 
                 inner join metadata_entry AS me_s on sme_s.metadata_id = me_s.id 
                 inner join sample_metadata_entry AS sme_c on sme_s.sample_id = sme_c.sample_id
                 inner join metadata_field AS mf_c on sme_c.metadata_KEY = mf_c.id
                 inner join metadata_entry AS me_c on sme_c.metadata_id = me_c.id 
                 inner join mlst_ecoli on me_c.value = mlst_ecoli.sample_code
               where me_s.value = '$serotype' and left(mf_s.label,9) = 'Antigen_O' and mf_c.label = 'Sample_code' and mlst_ecoli.permille_loci>799";
	}
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();

    open my $of, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
    while(my @row = $sth->fetchrow_array()){
      print $of "$row[0]\n"
    }       
    close $of;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
    return 0;
}