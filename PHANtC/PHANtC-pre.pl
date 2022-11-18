#!/usr/bin/env perl
## A wrapper script to obtain MLST profiles from IRIDA
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use DBI;
use Config::Simple;

# Parse arguments
my ($input,
    $species,
    $output) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
if ($species eq "Shiga toxin-producing Escherichia coli") { $species = "Escherichia coli"; };
my $idFastqs = getIdFiles($input);
createAllelesFile($idFastqs);
substituteCodesByNames($idFastqs);
adjustHeader();
move("cgMLST_header.tmp",$output);
exit(0);

# Obtain idFiles
sub getIdFiles {
    my ($forward) = @_;
    open my $if, '<', $forward;
    chomp(my $inIds = <$if>);
    close $if;
    return $inIds;
}

# Obtain allele profiles from db and write to temp file
sub createAllelesFile {
    my ($files_id) = @_;
    my $sql = "select allele_strain from v_mlst_ecoli_files_id where files_id = -1";
    if ($species eq "Escherichia coli") {
        $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union select allele_strain from v_mlst_ecoli_files_id where files_id IN ($idFastqs)";
    } else {
        if ($species eq "Listeria monocytogenes") {
            $sql = "select allele_strain from mlst_listeria where sample_code='FILE' union select allele_strain from v_mlst_listeria_files_id where files_id IN ($idFastqs)";
        }
    }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[0]\n";
    }       
    close $if;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
}

# Obtain samples names and codes from db and update the temp file
sub substituteCodesByNames {
    my ($files_id) = @_;
    my $sql = "select sampleName, sampleCode from v_sample_name_code where files_id IN ($files_id)";
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "code2name.tmp" or die "Cannot open code2name.tmp: $!";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[0]\t$row[1]\n";
    }       
    close $if;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
    # substitute the sample_code by its corresponding sample_name
    my $cmd = q( awk -F "\t" 'BEGIN {OFS = FS} NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' code2name.tmp cgMLST.tmp > cgMLST.tmp2 );
    system($cmd);
    move("cgMLST.tmp","cgMLST.tmp.OLD");
    move("cgMLST.tmp2","cgMLST.tmp");
}

# write the header in the EFSA way
sub adjustHeader {
    system("head -n 1 cgMLST.tmp > cgMLST_1.tmp");
    my $cmd = "echo ERROR";
    if ($species eq "Escherichia coli") {
        $cmd = q( awk 'BEGIN{FS="\t"; OFS=FS}{ for(i=2;i<=NF;i++){ $i = "INNUENDO_wgMLST-"$i".fasta" }}1' cgMLST_1.tmp > cgMLST_header.tmp );
    } else {
        if ($species eq "Listeria monocytogenes") {
            $cmd = q( awk 'BEGIN{FS="\t"; OFS=FS}{ for(i=2;i<=NF;i++){ $i = "Pasteur_cgMLST-"$i".fasta" }}1' cgMLST_1.tmp > cgMLST_header.tmp );
        }
    }	
    system($cmd);
    system("tail -n 1 cgMLST.tmp >> cgMLST_header.tmp");
}