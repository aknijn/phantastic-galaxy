#!/usr/bin/env perl
## A wrapper script to call PHANtA.py and collect its output
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use JSON::XS;

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
	# TRIMMING
    system("fastp --thread 4 -i $fastq1 -o fastq_1.fastq.gz -I $fastq2 -O fastq_2.fastq.gz -f 3 -t 3 -F 3 -T 3 -l 55 --cut_front_window_size 5 --cut_front_mean_quality 20 --cut_tail_window_size 5 --cut_tail_mean_quality 20");
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
    # JSON
    my $coverage = getCoverage();
	my ($read_mean_length, $q30_rate, $total_bases) = getFastp();
    open(my $fh, '>', $json) or die "Could not open file '$json' $!";
    print $fh "{\"coverage\": \"$coverage\",";
    print $fh "\"read_mean_length\": \"$read_mean_length\",";
    print $fh "\"q30_rate\": \"$q30_rate\",";
    print $fh "\"total_bases\": \"$total_bases\"}";
    close($fh);
    return 0;
}

sub getCoverage{
    my @trueCoverage_files = glob "output_dir/*/trueCoverage_report.txt";
	my $trueCoverage_file = @trueCoverage_files[0];
    open(my $fh, '<', $trueCoverage_file) or die "Could not open file '$trueCoverage_file' $!";
    my $lastline;
    $lastline = $_, while (<$fh>);
    chomp $lastline;
    close($fh);
    return $lastline;
}

sub getFastp{
    my $fname = 'fastp.json';
    my $txt   = do {                             
        local $/;                              
        open my $fh, "<", $fname or die $!;
        <$fh>;                                 
    };
    my $json = from_json($txt);
    my $read_mean_length = "$json->{summary}->{after_filtering}->{read1_mean_length}";
    my $q30_rate = "$json->{summary}->{after_filtering}->{q30_rate}";
    my $total_bases = "$json->{summary}->{after_filtering}->{total_bases}";
	return $read_mean_length, $q30_rate, $total_bases;
}

# Run QUAST
sub runQUAST {
    $genomeSize = $genomeSize * 1000000;
    system("quast --threads 4 -o outputdir --est-ref-size $genomeSize --min-contig 500 -l  \"$input_id\" --contig-thresholds 0,1000 $polished_fasta");
    system("mv outputdir/report.tsv $quast");
    return 0;
}
