#!/usr/bin/env perl
## A wrapper script to call MentaLiST
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
    $input2,
    $input1_name,
	$phantclm_allele,
	$phantclm_json,
    $phantclm_tree,
    $phantclm_dm) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $sampleGenesMapped = "0";
my $sampleGenesNovel = "0";
my $sampleGenesMultiple = "0";
my $sampleGenesPartial = "0";
my $sampleGenesNotMapped = "0";
runMentaLiST();
exit(0);

# Run runMentaLiST
sub runMentaLiST {
    my $mlst_schema = "/afs/galaxy/tool-data/mentalist_databases/listeria_monocytogenes_moura_k31_2018-10-08/listeria_monocytogenes_moura_k31_m023_2018-10-08.jld";
	open my $ipf, '>', 'inputfile' or die "Cannot open inputfile: $!";
    print $ipf "$input1_name\t$input1\n";
    if ($input2 ne "NULL") {
      print $ipf "$input1_name\t$input2\n";
    }
	close $ipf;
	system("mentalist call --output_votes -o mentalist_out --db $mlst_schema -i inputfile");
    # add mentalist_out to database
    my $permille_loci = collectOutput();
    createAllelesFile();
    # get statistics
    open my $f1, '<', "cgMLST.tmp"; 
    my $f1Line = <$f1>; 
    close $f1;
    my @f1columns = split('\t', $f1Line);
    my $cgLociNumber = scalar @f1columns - 3;
	$sampleGenesMapped = `grep -c "Called" mentalist_out.coverage.txt`;
	chomp $sampleGenesMapped;
	$sampleGenesNovel = `grep -c "Novel" mentalist_out.coverage.txt`;
	chomp $sampleGenesNovel;
	$sampleGenesMultiple = `grep -c "Multiple" mentalist_out.coverage.txt`;
	chomp $sampleGenesMultiple;
	$sampleGenesPartial = `grep -c "Partially" mentalist_out.coverage.txt`;
	chomp $sampleGenesPartial;
	$sampleGenesNotMapped = `grep -c "Not present" mentalist_out.coverage.txt`;
	chomp $sampleGenesNotMapped;
    open my $of, '>', $phantclm_json or die "Cannot open json: $!";
    print $of "{\"core_genome_schema_size\": $cgLociNumber, \"sample_genes_mapped\": $sampleGenesMapped}";
    close $of;
    open my $of2, '>', $phantclm_allele or die "Cannot open allele: $!";
    print $of2 "Genome\tEXC\tLNF\tINF\tMULT\tPART\n";
    print $of2 "$input1_name\t$sampleGenesMapped\t$sampleGenesNotMapped\t$sampleGenesNovel\t$sampleGenesMultiple\t$sampleGenesPartial";
    close $of2;
    # calc distance matrix from database
    my $result = system("python $scriptdir/scripts/mlst_hash_stretch_distance.py -i cgMLST.tmp -o $phantclm_dm");
    # calc tree from distance matrix
    system("$scriptdir/scripts/mentalist_tree $phantclm_dm > $phantclm_tree");
    return 0;
}

# Collect output to database
sub collectOutput{
    my $alleles = "mentalist_out";
    open(my $if_in, $alleles) or die "Could not read from mentalist_out, program halting.";
    <$if_in>;
    my $allele_line = <$if_in>;
    chomp $allele_line;
    close $if_in;
    my @loci = split(/\t/, $allele_line);
	my $sample_code = $loci[0];
	my $numloci = scalar @loci - 3;
	my $mappedloci = `grep -c "Called" mentalist_out.coverage.txt`;
	my $permille_loci = int(($mappedloci/$numloci)*1000 + 0.5);

    my $sql_insert = "insert into mlst_listeria (sample_code) values (?)";
    my $sql_update = "update mlst_listeria set permille_loci=?, allele_strain=? where sample_code=?";
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
    return $permille_loci;
}

# Obtain alleles from db with allele_coverage>79,9% and write to temp file
sub createAllelesFile {
    my $sql = "select mlst_listeria.allele_strain from mlst_listeria where permille_loci>799";
    # connect to MySQL database
    my %attr2 = ( PrintError=>0, RaiseError=>1);
    my $dbh2 = DBI->connect($dsn,$user,$pwd,\%attr2);
    my $sth2 = $dbh2->prepare($sql);
    $sth2->execute();

    open my $if, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
    while(my @row = $sth2->fetchrow_array()){
      print $if "$row[0]\n"
    }       
    close $if;

    $sth2->finish();
    # disconnect from the MySQL database
    $dbh2->disconnect();
    return 0;
}