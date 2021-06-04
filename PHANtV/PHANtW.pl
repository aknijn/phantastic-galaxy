#!/usr/bin/env perl
## A wrapper script to call PopPUNK.py
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;

# Parse arguments
my ($input_file,
    $phantw_cl,
    $phantw_tree) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
prepareEnvironment();
my $result=runPopPUNK();
move("output/output_core_NJ.nwk", $phantw_tree);
system("sed 's-input_dir/--g' distances.txt > $phantw_cl");

exit($result);

# Run PopPUNK
sub runPopPUNK {
    my $result_db = system("poppunk --create-db --r-files inputfiles.txt --threads 4 --k-step 2 --min-k 9 --plot-fit 0 --overwrite --output dbdir");
    my $result_fit = system("poppunk --fit-model --threads 4 --output output --full-db --K 2 --microreact --ref-db dbdir --distances dbdir/dbdir.dists");
    system("python $scriptdir/scripts/extract_distances.py  --distances dbdir/dbdir.dists --output distances.txt");
    return $result_db+$result_fit;
}

# Run prepareEnvironment, create the directory $indir with symlinks to the files listed in $input_file
sub prepareEnvironment {
    my $indir = "input_dir";
    my $output_file = "inputfiles.txt";
    my $line = "";
    mkdir($indir);
    open(my $if, $input_file) or die "Could not read from input_file, program halting.";
    open(my $of, ">", $output_file) or die "Error writing to output_file, program halting.";
    while ($line = <$if>) {
        chomp $line;
        my ($idSample, $fastafile) = split(/\t/, $line);
        if ($fastafile ne "NULL") {
            symlink($fastafile, $indir . "/" . $idSample);
            print $of $indir . "/" . $idSample . "\n";
        }
    }
    close $if;
    close $of;
    return 0;
}


