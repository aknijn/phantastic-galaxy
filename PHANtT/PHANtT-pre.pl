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
my $iridadir = $cfg->param('fs.output_path');
my $idFastqs = getIdFiles($input);
getSampleNames($idFastqs);
exit(0);

sub getSampleNames{
    my $idFastqs = $_[0];
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sql;

    $sql = "select sampleName, file_path from v_sample_name_virulotypes where files_id IN ($idFastqs)";

    my $sth = $dbh->prepare($sql);
    $sth->execute();
  
    open my $if, '>', $output or die "Error writing to output, program halting.";
    while (my @row = $sth->fetchrow_array) { 
        print $if "$row[0]\t/$iridadir/$row[1]\n";
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
