#!/usr/bin/env perl
## A wrapper script to obtain sample names from IRIDA
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
    $output) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $idFastqs = getIdFiles($input);
getSampleNames($idFastqs);
exit(0);

sub getSampleNames{
    my $idFastqs = $_[0];
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sql;

    $sql = "select
             sample.sampleName, analysis_output_file.file_path
           from
             sequence_file_pair_files
             inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id
             inner join sample on sample_sequencingobject.sample_id = sample.id
             inner join analysis_submission_sequencing_object on sequence_file_pair_files.pair_id = analysis_submission_sequencing_object.sequencing_object_id
             inner join analysis_submission on analysis_submission_sequencing_object.analysis_submission_id = analysis_submission.id
             inner join analysis_output_file_map on analysis_submission.analysis_id = analysis_output_file_map.analysis_id
             inner join analysis_output_file on analysis_output_file_map.analysisOutputFilesMap_id = analysis_output_file.id
           where analysis_output_file_map.analysis_output_file_key = 'virulotypes' and sequence_file_pair_files.files_id IN ($idFastqs)";

    my $sth = $dbh->prepare($sql);
    $sth->execute();
  
    open my $if, '>', $output or die "Error writing to output, program halting.";
    while (my @row = $sth->fetchrow_array) { 
        print $if "$row[0]\t/afs/irida/data/output/$row[1]\n";
    }       
    close $if;
       
    $sth->finish();
     # disconnect from the MySQL database
    $dbh->disconnect();
}

# Obtain idFiles from file paths
sub getIdFiles {
    my ($forward) = @_;
    open my $if, '<', $forward;
    chomp(my $inIds = <$if>);
    close $if;
    return $inIds;
}
