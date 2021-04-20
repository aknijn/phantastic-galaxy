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
    my $prepath = "/afs/irida21/data/output/";
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sql;
    if ($useNames eq "true") 
    {
        $sql = "select sampleName, file_path from v_sample_name_code_fasta where files_id IN ($idFastqs)";
    } else {
        $sql = "select sampleCode, file_path from v_sample_name_code_fasta where files_id IN ($idFastqs)";
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