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
    $seqtype,
    $python) = @ARGV;

# Run program
if (index($input1, '.fasta') >= 0) { $input1 = "NULL" }
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $result = runLisSero();
runMLST();
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

# Run MLST
sub runMLST {
    system("mlst --legacy --scheme listeria_2 " . $inputf . " | cut -f3,4,5,6,7,8,9,10 > mlstsevenloci");
    copy("mlstsevenloci",$seqtype) or die "Could not copy mlstsevenloci, program halting.";
    return 0;
}

# Run Virulotyper
sub runVirulotyper {
    if ($input1 ne "NULL") {
      if ($input2 ne "NULL") {
        system("perl $scriptdir/scripts/patho_typing.pl 'python $scriptdir/scripts/patho_typing.py -s Listeria monocytogenes -f $input1 $input2 -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'");
        system("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > $virulotypes");
      } else {
        system("perl $scriptdir/scripts/patho_typing.pl 'python $scriptdir/scripts/patho_typing.py -s Listeria monocytogenes -f $input1 -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'");
        system("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > $virulotypes");
      }
    } else {
      system("touch $virulotypes");
    }
    return 0;
}

# Run AMRgenes
sub runAMRgenes {
    if ($input1 ne "NULL") {
      system("abricate --db ncbi $inputf > $amrgenes");
    } else {
      system("touch $amrgenes");
    }
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
    # MLST
    my $sequence_qc_result = 1;
    $input_file = 'mlstsevenloci';
    open($if, $input_file) or die "Could not read from mlstsevenloci, program halting.";
    <$if>;
    $line = <$if>;
    chomp $line;
    my @mlst_st = split(/\t/, $line);
    close $if;
    if ((index($line, '?') != -1) or (index($line, '-') != -1)) { $sequence_qc_result = 0; }
    my $qc_status = "Failed";
    my $qc_messages = "Accepted for outbreak investigation.";
    if ($sequence_qc_result) {
        $qc_status = "Passed";
        $qc_messages = "Passed.";
    }
    my @CCandLineage = getCCandLineage($mlst_st[0]);
    # write json
    open(my $of, ">", $output_file) or die "Could not read from output_tab, program halting.";
    print $of "\{\"information_name\": \"" . $input_name . "\", \"qc_status\": \"" . $qc_status . "\", \"qc_messages\": \"" . $qc_messages . "\", \"serotype_serogroup\": \"" . $serotype . "\", \"serotype_amplicons\": \"" . $amplicons . "\", \"mlst_ST\": \"ST" . $mlst_st[0] .  "\", \"mlst_CC\": \"" . $CCandLineage[0] . "\", \"mlst_lineage\": \"" . $CCandLineage[1] . "\", \"region\": \"" . $region . "\", \"year\": \"" . $year . "\"\}";
    close $of;
    return 0;
}

sub getCCandLineage {
    my $line;
    my ($st) = @_;
    my @result = ("-", "-");
    my $mlst_profile = "$scriptdir/data/lmonocytogenes.txt";
    open my $if, '<', $mlst_profile;
    <$if>;
    while ($line = <$if>) {
        chomp $line;
        my @profile = split(/\t/, $line);
        if ($profile[0] eq $st) { @result = ($profile[8], $profile[9]); }
    }
    return @result;
}
