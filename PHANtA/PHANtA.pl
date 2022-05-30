#!/usr/bin/env perl
## A wrapper script to call PHANtA.py and collect its output
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use IO::Compress::Gzip qw(gzip $GzipError) ;

# Parse arguments
my ($fastq1,
    $fastq2,
    $input_id,
    $genomeSize,
    $polished_fasta,
    $json,
    $quast,
    $python) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
runPHANtA();
collectOutput();
runQUAST();
exit 0;

# Run PHANtA
sub runPHANtA {
    gzip $fastq1 => "fastq_1.fastq.gz" or die "gzip failed: $GzipError\n";
    gzip $fastq2 => "fastq_2.fastq.gz" or die "gzip failed: $GzipError\n";
    my $rematchdir = "$scriptdir/../ReMatCh";
    my $mode = 0744;
    chmod $mode, "$rematchdir/rematch.py";
    my $newpath = "PATH=$ENV{PATH}:$rematchdir";
    `$newpath; $python`;
    return 0;
}

# Collect output
sub collectOutput{
    # FASTA
    my $final_file = 'output_dir/fastq/final_assembly.txt';
    open(my $if, $final_file) or die "Could not read from final_assembly.txt, program halting.";
    my $line = <$if>;
    chomp $line;
    move($line, $polished_fasta);
    close $if;
    # COVERAGE
    my @trueCoverage_files = glob "output_dir/*/trueCoverage_report.txt";
    foreach my $trueCoverage_file (@trueCoverage_files)
    {
      open(my $fh, '<', $trueCoverage_file) or die "Could not open file '$trueCoverage_file' $!";
      my $lastline;
      $lastline = $_, while (<$fh>);
      chomp $lastline;
      close($fh);
      open($fh, '>', $json) or die "Could not open file '$json' $!";
      print $fh "{\"coverage\": \"$lastline\"}";
      close($fh);
    }
    return 0;
}

# Run QUAST
sub runQUAST {
    $genomeSize = $genomeSize * 1000000;
    system("quast --threads 4 -o outputdir --est-ref-size $genomeSize --min-contig 500 -l  \"$input_id\" --contig-thresholds 0,1000 $polished_fasta");
    system("mv outputdir/report.tsv $quast");
    return 0;
}




