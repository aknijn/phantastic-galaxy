#!/usr/bin/env perl
## A wrapper script to call LisSero.py
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;

# Parse arguments
my ($input1,
    $input2,
    $inputf,
    $input_name,
    $region,
    $year,
    $output_file,
    $virulotypes,
    $amrgenes,
    $python) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $result = runLisSero();
runMentaLiST();
runVirulotyper();
runAMRgenes();
collectOutput();
exit($result);

# Run LisSero
sub runLisSero {
    $python = $python . " " . $inputf . " > output_tab";
    print "Running $python \n";
    my $result = system("$python");
    return $result;
}

# Run runMentaLiST
sub runMentaLiST {
    my $mlst_schema = "$scriptdir/data/listeria_monocytogenes_pubmlst.jld";
    if ($input2 ne "NULL") {
      system("mentalist call --output_votes -o mentalist_out --db $mlst_schema -1 $input1 -2 $input2");
    } else {
      system("mentalist call --output_votes -o mentalist_out --db $mlst_schema -1 $input1");
    }
    return 0;
}

# Run Virulotyper
sub runVirulotyper {
    if ($input2 ne "NULL") {
      system("perl $scriptdir/scripts/patho_typing.pl 'python $scriptdir/scripts/patho_typing.py -s Listeria monocytogenes -f $input1 $input2 -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'");
      system("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > $virulotypes");
    } else {
      system("perl $scriptdir/scripts/patho_typing.pl 'python $scriptdir/scripts/patho_typing.py -s Listeria monocytogenes -f $input1 -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'");
      system("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > $virulotypes");
    }
    return 0;
}

# Run AMRgenes
sub runAMRgenes {
    system("abricate --db ncbi $inputf > $amrgenes");
    return 0;
}

# Collect output
sub collectOutput{
    # LisSero
    my $input_file = 'output_tab';
    open(my $if, $input_file) or die "Could not read from output_tab, program halting.";
    <$if>;
    my $line = <$if>;
    chomp $line;
    my ($strain, $serotype, $amplicons) = split(/\t/, $line);
    close $if;
    # MentaLiST
    $input_file = 'mentalist_out.byvote';
    open($if, $input_file) or die "Could not read from mentalist_out, program halting.";
    <$if>;
    $line = <$if>;
    chomp $line;
    my @mlst_st = split(/\t/, $line);
    close $if;
    # Quality Control
    my $sequence_qc_result = 1;
    my $sequence_qc_cov = 0;
    $input_file = 'mentalist_out.coverage.txt';
    open $if, '<', $input_file;
    <$if>;
    while ($line = <$if>) {
        chomp $line;
        my @sequence_qc = split(/\t/, $line);
        if ($sequence_qc[2] < 1) { $sequence_qc_result = 0; }
        $sequence_qc_cov = $sequence_qc_cov + $sequence_qc[3];
    }
    if ($sequence_qc_cov < 210) { $sequence_qc_result = 0; }
    close $if;
    my $qc_status = "Failed";
    my $qc_messages = "Accepted for outbreak investigation.";
    if ($sequence_qc_result) {
        $qc_status = "Passed";
        $qc_messages = "Passed.";
    }
    my $lineage = getLineage($mlst_st[8]);
    # write json
    open(my $of, ">", $output_file) or die "Could not read from output_tab, program halting.";
    print $of "\{\"information_name\": \"" . $input_name . "\", \"qc_status\": \"" . $qc_status . "\", \"qc_messages\": \"" . $qc_messages . "\", \"serotype_serogroup\": \"" . $serotype . "\", \"serotype_amplicons\": \"" . $amplicons . "\", \"mlst_ST\": \"ST" . $mlst_st[8] .  "\", \"mlst_CC\": \"" . $mlst_st[9] . "\", \"mlst_lineage\": \"" . $lineage . "\", \"region\": \"" . $region . "\", \"year\": \"" . $year . "\"\}";
    close $of;
    return 0;
}

sub getLineage {
    my $line;
    my ($st) = @_;
    my $result = "-";
    my $mlst_profile = "$scriptdir/data/listeria_monocytogenes_pubmlst.jld.profile";
    open my $if, '<', $mlst_profile;
    <$if>;
    while ($line = <$if>) {
        chomp $line;
        my @profile = split(/\t/, $line);
        if ($profile[0] eq $st) { $result = $profile[9]; }
    }
    return $result;
}
