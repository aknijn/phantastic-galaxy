#!/usr/bin/env perl
## A wrapper script to collect patho_typing.py output
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;

# Parse arguments
my ($python) = @ARGV;

# Run program
runPathoTyping();
collectOutput();
exit 0;

# Run patho_typing
sub runPathoTyping {
    my $abs_path = Cwd::abs_path($PROGRAM_NAME);
    my $scriptdir = dirname($abs_path);
    my $rematchdir = "$scriptdir/../../ReMatCh";
    my $newpath = "PATH=$ENV{PATH}:$rematchdir";
    `$newpath; $python`;
     return 0;
}

sub collectOutput{
    my $patho_type = "";
    open(my $fh, '<', 'output_dir/rematch/rematchModule_report.txt') || die "Could not open file 'output_dir/rematch/rematchModule_report.txt' $!";
    my @rematch_lines = <$fh>;
    close($fh);
    my @rematch_total_table = "";
    foreach(@rematch_lines) {
      if ($_ =~ m/\t/) { # Only table lines
        my @elems = split('\t', $_);
        if ($_ =~ m/#/) { # First line
          s/_/ /g for @elems;
          @rematch_total_table = join("\t", @elems[0,1,2,5]);
        }
        else {
          push(@rematch_total_table,join("\t", @elems[0,1,2,5]));
        }
      }
    }
    open($fh, '>', 'pathotyper_rep_tot_tab') || die "Could not open file 'pathotyper_rep_tot_tab' $!";
    print $fh @rematch_total_table;
    close($fh);

    return 0;
}
