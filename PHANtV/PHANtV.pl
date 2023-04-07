#!/usr/bin/env perl
## A wrapper script to call allele distance matrix and tree
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use DBI;
use Config::Simple;

# Parse arguments
my ($input_files,
    $species,
    $useNames,
    $phantv_am,
    $phantv_dm,
    $phantv_tree) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
if ($species eq "Shiga toxin-producing Escherichia coli") { $species = "Escherichia coli"; };
runAlleleObserver();
exit(0);

# Get allele file with clean columns and calculate the distance matrix and the phylogenetic NJ tree
sub runAlleleObserver {
    # get idFiles from filepaths
    open my $if, '<', $input_files;
    chomp(my @inputs = <$if>);
    close $if;
    my $files_id = getIdFilesString(@inputs);
    createAllelesFile($files_id);
    if ($useNames eq "true") { substituteCodesByNames($files_id); }
    # copy filtered allele profiles in allele matrix file (columns with 0 are removed)
    my $cmd = q( awk -F"\t" '{for(i=2;i<=NF;i++){if ($i == 0){print i}}}' cgMLST.tmp | sort | uniq > cgMLST.nocols );
    system($cmd);
    open my $tf, '<', 'cgMLST.nocols';
    chomp(my @noCols = <$tf>);
    close $tf;
    my $noColsString = join ',', @noCols;
	if ($noColsString eq "") {
        copy("cgMLST.tmp",$phantv_am);
    } else {
        my $cmd2 = "cut --complement -f" . $noColsString . " cgMLST.tmp > " . $phantv_am;
        system($cmd2);
    }
    # calculate distance matrix from allele profiles
    my $result = system("python $scriptdir/scripts/mlst_hash_stretch_distance.py -i cgMLST.tmp -o $phantv_dm");
    # calculate tree from distance matrix
    system("python $scriptdir/scripts/mentalist_tree.py $phantv_dm > $phantv_tree");
    return 0;
}

# Obtain idFile from file path
sub getIdFile {
    my ($inpath) = @_;
    my(@dirs) = split m%/%, $inpath;
    return $dirs[5];
}

# Obtain a string with the idFiles
sub getIdFilesString {
    my @inFiles = @_;
    my @inIdFiles;
    foreach (@inFiles) {
        push @inIdFiles, getIdFile($_);
    }
    my $idFilesString = join ',', @inIdFiles;
    return $idFilesString;
}

# Obtain allele profiles from db and write to temp file
sub createAllelesFile {
    my ($files_id) = @_;
    my $sql;

    if (index($files_id, ',') == -1) {
        # only one file in input, return allele strains of samples with same serotype
        if ($species eq "Escherichia coli") {
            $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
                    select allele_strain from v_mlst_ecoli_files_id_1 where files_id = $files_id";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select mlst_listeria.allele_strain from mlst_listeria";
            }
        }
    } else {
        if ($species eq "Escherichia coli") {
            $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
                    select allele_strain from v_mlst_ecoli_files_id where files_id IN ($files_id)";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select allele_strain from mlst_listeria where sample_code='FILE' union
                        select allele_strain from v_mlst_listeria_files_id where files_id IN ($files_id)";
            }
        }
    }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
	# print $if "FILE\t";
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
    my $sql;
    if (index($files_id, ',') == -1) {
        # only one file in input, return allele strains of samples with same serotype
        if ($species eq "Escherichia coli") {
            $sql = "select sampleName, sampleCode from v_sample_name_code_ecoli_1 where files_id IN ($files_id)";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select sampleName, sampleCode from v_sample_name_code_listeria_1";
            }
        }
    } else { $sql = "select sampleName, sampleCode from v_sample_name_code where files_id IN ($files_id)"; }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "code2name.tmp" or die "Cannot open code2name.tmp: $!";
	print $if "FILE\tFILE\n";
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
