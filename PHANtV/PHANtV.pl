#!/usr/bin/env perl
## A wrapper script to call MentaLiST distance matrix and tree
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
runMentaLiST();
exit(0);

# Run runMentaLiST
sub runMentaLiST {
    # get idFiles from filepaths
    open my $if, '<', $input_files;
    chomp(my @inputs = <$if>);
    close $if;
    my $files_id = getIdFilesString(@inputs);
    createAllelesFile($files_id);
    if ($useNames eq "true") { substituteCodesByNames($files_id); }
    # calc distance matrix from temp database
    #system("$scriptdir/scripts/mentalist_distance cgMLST.tmp > $phantv_dm");
    my $result = system("python $scriptdir/scripts/mlst_hash_stretch_distance.py -i cgMLST.tmp -o $phantv_dm");
    # calc tree from distance matrix
    system("$scriptdir/scripts/mentalist_tree $phantv_dm > $phantv_tree");
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

# Obtain alleles from db and write to temp file
sub createAllelesFile {
    my ($files_id) = @_;
    my $sql;

    if (index($files_id, ',') == -1) {
        # only one file in input, return allele strains of samples with same serotype
        if ($species eq "Escherichia coli") {
            $sql = "select mlst_ecoli.allele_strain from mlst_ecoli where sample_code='FILE' union
                    select mlst_ecoli.allele_strain from sequence_file_pair_files AS sfpf_s
                      inner join sample_sequencingobject AS sso_s on sfpf_s.pair_id = sso_s.sequencingobject_id
                      inner join sample_metadata_entry AS sme_s on sso_s.sample_id = sme_s.sample_id
                      inner join metadata_field AS mf_s on sme_s.metadata_KEY = mf_s.id
                      inner join metadata_entry AS me_s on sme_s.metadata_id = me_s.id 
                      inner join metadata_entry AS me_m on me_s.value = me_m.value 
                      inner join sample_metadata_entry AS sme_m on me_m.id = sme_m.metadata_id
                      inner join sample_metadata_entry AS sme_c on sme_m.sample_id = sme_c.sample_id
                      inner join metadata_field AS mf_c on sme_c.metadata_KEY = mf_c.id
                      inner join metadata_entry AS me_c on sme_c.metadata_id = me_c.id 
                      inner join mlst_ecoli on me_c.value = mlst_ecoli.sample_code
                    where left(mf_s.label,9) = 'Antigen_O' and mf_c.label = 'Sample_code' and sfpf_s.files_id = $files_id";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select mlst_listeria.allele_strain from mlst_listeria";
            }
        }
    } else {
        if ($species eq "Escherichia coli") {
            $sql = "select mlst_ecoli.allele_strain from mlst_ecoli where sample_code='FILE' union
                    select mlst_ecoli.allele_strain from sequence_file_pair_files
                      inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
                      inner join sample_metadata_entry on sample_sequencingobject.sample_id = sample_metadata_entry.sample_id
                      inner join metadata_field on sample_metadata_entry.metadata_KEY = metadata_field.id
                      inner join metadata_entry on sample_metadata_entry.metadata_id = metadata_entry.id
                      inner join mlst_ecoli on metadata_entry.value = mlst_ecoli.sample_code
                    where left(metadata_field.label,11) = 'Sample_code' and sequence_file_pair_files.files_id IN ($files_id)";} else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select mlst_listeria.allele_strain from mlst_listeria where sample_code='FILE' union
                        select mlst_listeria.allele_strain from sequence_file_pair_files
                          inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
                          inner join sample_metadata_entry on sample_sequencingobject.sample_id = sample_metadata_entry.sample_id
                          inner join metadata_field on sample_metadata_entry.metadata_KEY = metadata_field.id
                          inner join metadata_entry on sample_metadata_entry.metadata_id = metadata_entry.id
                          inner join mlst_listeria on metadata_entry.value = mlst_listeria.sample_code
                        where metadata_field.label = 'Sample_code' and sequence_file_pair_files.files_id IN ($files_id)";}
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
    my $sql = "select sample.sampleName, metadata_entry.value from sequence_file_pair_files
                  inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
                  inner join sample_metadata_entry on sample_sequencingobject.sample_id = sample_metadata_entry.sample_id
                  inner join metadata_field on sample_metadata_entry.metadata_KEY = metadata_field.id
                  inner join metadata_entry on sample_metadata_entry.metadata_id = metadata_entry.id
                  inner join sample on sample_sequencingobject.sample_id = sample.id
                where left(metadata_field.label,11) = 'Sample_code' and sequence_file_pair_files.files_id IN ($files_id)";

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
    # substitute the sample_code by its corrisponding sample_name
    my $cmd = q( awk -F "\t" 'BEGIN {OFS = FS} NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' code2name.tmp cgMLST.tmp > cgMLST.tmp2 );
    system($cmd);
    system("mv cgMLST.tmp cgMLST.tmp.OLD");
    system("mv cgMLST.tmp2 cgMLST.tmp");
}
