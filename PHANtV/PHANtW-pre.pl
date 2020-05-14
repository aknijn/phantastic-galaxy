#!/usr/bin/env perl
## A wrapper script to call PopPUNK.py
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
    $useNames,
    $output_file) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $idFastqs = getIdFiles($input1);
getFastaPaths($idFastqs);
exit(0);

sub getFastaPaths{
    my $idFastqs = $_[0];
    my $prepath = "/afs/irida/data/output/";
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sql;
    if ($useNames eq "true") 
	{
        $sql = "select
                 sample.sampleName, analysis_output_file.file_path
               from
                 sequence_file_pair_files
                 inner join analysis_submission_sequencing_object on sequence_file_pair_files.pair_id = analysis_submission_sequencing_object.sequencing_object_id
                 inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
                 inner join analysis_submission on analysis_submission_sequencing_object.analysis_submission_id = analysis_submission.id
                 inner join analysis_output_file_map on analysis_submission.analysis_id = analysis_output_file_map.analysis_id
                 inner join analysis_output_file on analysis_output_file_map.analysisOutputFilesMap_id = analysis_output_file.id
                 inner join sample on sample_sequencingobject.sample_id = sample.id
               where analysis_output_file_map.analysis_output_file_key = 'phantastic_contigs' and sequence_file_pair_files.files_id IN ($idFastqs)";
    } else {
        $sql = "select
                 metadata_entry.value, analysis_output_file.file_path
               from
                 sequence_file_pair_files
                 inner join analysis_submission_sequencing_object on sequence_file_pair_files.pair_id = analysis_submission_sequencing_object.sequencing_object_id
                 inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
                 inner join analysis_submission on analysis_submission_sequencing_object.analysis_submission_id = analysis_submission.id
                 inner join analysis_output_file_map on analysis_submission.analysis_id = analysis_output_file_map.analysis_id
                 inner join analysis_output_file on analysis_output_file_map.analysisOutputFilesMap_id = analysis_output_file.id
                 inner join sample_metadata_entry on sample_sequencingobject.sample_id = sample_metadata_entry.sample_id
                 inner join metadata_field on sample_metadata_entry.metadata_KEY = metadata_field.id
                 inner join metadata_entry on sample_metadata_entry.metadata_id = metadata_entry.id
               where analysis_output_file_map.analysis_output_file_key = 'phantastic_contigs' and left(metadata_field.label,11) = 'Sample_code' and sequence_file_pair_files.files_id IN ($idFastqs)";
	}
    my $sth = $dbh->prepare($sql);
    $sth->execute();
  
    open my $if, '>', $output_file or die "Error writing to output_file, program halting.";
    while (my @row = $sth->fetchrow_array) { 
        print $if "$row[0]\t$prepath$row[1]\n";
    }       
    close $if;
       
    $sth->finish();
     # disconnect from the MySQL database
    $dbh->disconnect();
}

# Obtain idFiles from file paths
sub getIdFiles {
    my ($inpath) = @_;
    open my $if, '<', $inpath;
    chomp(my @inFiles = <$if>);
    close $if;
    my @inIdFiles;
    foreach (@inFiles) {
        push @inIdFiles, getIdFile($_);
    }
    my $files_id = join ',', @inIdFiles;
    return $files_id;
}

# Obtain idFile from file path
sub getIdFile {
    my ($inpath) = @_;
    my(@dirs) = split m%/%, $inpath;
    return $dirs[5];
}