#!/usr/bin/env perl
## A wrapper script to call grapetree tree with metadata
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use File::Path;
use DBI;
use Config::Simple;

# Parse arguments
my ($input_files,
    $species,
    $phantg_am,
    $phantg_metadata,
    $phantg_tree,
    $phantg_grapetree) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
if ($species eq "Shiga toxin-producing Escherichia coli") { $species = "Escherichia coli"; };
my $files_id = getIdFilesString($input_files);
runGrapeTree($files_id);
createMetadataFile($files_id);
createGrapeTreeLink();
exit(0);

# Get allele file with clean columns and calculate the phylogenetic Minimum Spanning Tree using GrapeTree
sub runGrapeTree {
    my ($files_id) = @_;
    createAllelesFile($files_id);
    substituteCodesByNames($files_id);
    # copy filtered allele profiles in allele matrix file (columns with INF, LNF, ecc. are removed)
    my $cmd = q( awk -F"\t" '{(NR>1)}{for(i=2;i<=NF;i++){if ($i ~ /^0$|[a-zA-Z+]+/){print i}}}' cgMLST.tmp | sort | uniq > cgMLST.nocols );
    system($cmd);
    open my $tf, '<', 'cgMLST.nocols';
    chomp(my @noCols = <$tf>);
    close $tf;
    my $noColsString = join ',', @noCols;
    if ($noColsString eq "") {
        copy("cgMLST.tmp",$phantg_am);
    } else {
        $cmd = "cut --complement -f" . $noColsString . " cgMLST.tmp > " . $phantg_am;
        system($cmd);
    }
    # calculate tree from allele profiles
    $cmd = q( ln -s $(dirname $(which grapetree))/../lib/python3.8/site-packages/grapetree/grapetree.py grapetree.py );
    system($cmd);
    my $result = system("python grapetree.py -p $phantg_am -m MSTreeV2 > $phantg_tree");
    return 0;
}

# Obtain idFile from file path
sub getIdFile {
    my ($inpath) = @_;
    my(@dirs) = split m%/%, $inpath;
    return $dirs[5];
}

# Obtain a string with the idFiles from filepaths
sub getIdFilesString {
    my ($inputFiles) = @_;
    open my $if, '<', $inputFiles;
    chomp(my @inFiles = <$if>);
    close $if;
    my @inIdFiles;
    foreach (@inFiles) {
        push @inIdFiles, getIdFile($_);
    }
    my $idFilesString = join ',', @inIdFiles;
    return $idFilesString;
}

# Obtain allele profiles from db and write to temp file
sub createAllelesFile {
    my ($files_id) = @_;
    my $sql;

    if (index($files_id, ',') == -1) {
        # only one file in input, return allele strains of samples with same serotype
        if ($species eq "Escherichia coli") {
            $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
                    select allele_strain from v_mlst_ecoli_files_id_1 where files_id = $files_id";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select mlst_listeria.allele_strain from mlst_listeria";
            }
        }
    } else {
        if ($species eq "Escherichia coli") {
            $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
                    select allele_strain from v_mlst_ecoli_files_id where files_id IN ($files_id)";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select allele_strain from mlst_listeria where sample_code='FILE' union
                        select allele_strain from v_mlst_listeria_files_id where files_id IN ($files_id)";
            }
        }
    }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "cgMLST.tmp" or die "Cannot open cgMLST.tmp: $!";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[0]\n";
    }       
    close $if;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
}

# Obtain samples names and codes from db and update the temp file
sub substituteCodesByNames {
    my ($files_id) = @_;
    my $sql;
    if (index($files_id, ',') == -1) {
        # only one file in input, return allele strains of samples with same serotype
        if ($species eq "Escherichia coli") {
            $sql = "select sampleName, sampleCode from v_sample_name_code_ecoli_1 where files_id IN ($files_id)";
        } else {
            if ($species eq "Listeria monocytogenes") {
                $sql = "select sampleName, sampleCode from v_sample_name_code_listeria_1";
            }
        }
    } else { $sql = "select sampleName, sampleCode from v_sample_name_code where files_id IN ($files_id)"; }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "code2name.tmp" or die "Cannot open code2name.tmp: $!";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[0]\t$row[1]\n";
    }       
    close $if;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
    # substitute the sample_code by its corresponding sample_name
    my $cmd = q( awk -F "\t" 'BEGIN {OFS = FS} NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' code2name.tmp cgMLST.tmp > cgMLST.tmp2 );
    system($cmd);
    move("cgMLST.tmp","cgMLST.tmp.OLD");
    move("cgMLST.tmp2","cgMLST.tmp");
}

# Obtain metadata from db and write to output file
sub createMetadataFile {
    my ($files_id) = @_;
    my $sql;

    if ($species eq "Escherichia coli") {
        $sql = "select * from v_grapetree_ecoli where files_id IN ($files_id)";
    } else {
        if ($species eq "Listeria monocytogenes") {
            $sql = "select * from v_grapetree_listeria where files_id IN ($files_id)";
        } else { $sql = "select * from v_grapetree_listeria where files_id = -1"; }
    }

    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "$phantg_metadata" or die "Cannot open $phantg_metadata: $!";
    if ($species eq "Escherichia coli") { print $if "Campione\tRegione\tAnno\tAntigen_O\tAntigen_H\tCluster_id\teae\tehxA\tMLST_ST\tQC\tstx1\tstx2\tstx_subtype\n"; } 
    else { if ($species eq "Listeria monocytogenes") { print $if "Campione\tRegione\tAmplicons\tAnno\tCluster_id\tMLST_CC\tMLST_Lineage\tQC\tSeroGroup\tSeroType\n"; } }
    while (my @row = $sth->fetchrow_array) { 
      if ($species eq "Escherichia coli") { print $if "$row[1]\t$row[10]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[11]\t$row[12]\t$row[13]\n"; } 
      else { if ($species eq "Listeria monocytogenes") { print $if "$row[1]\t$row[8]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[9]\t$row[10]\n"; } }
    }       
    close $if;

    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
}

sub createGrapeTreeLink {
    # copy $phantg_metadata & $phantg_tree to a shared path visible by GrapeTree on the server vizapp.iss.it
    my $grape_path = "/gfs/vizapp/grapetree/";
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year = 1900 + $year;
    my $stamp = sprintf("%04d%02d%02d%02d%02d%02d", $year, $mon+1, $mday, $hour, $min, $sec);
    my $grape_path_yearmm = $grape_path . substr($stamp,0,6);
    # create the path if it doesn't exist yet
    if ( !-d $grape_path_yearmm) { mkdir $grape_path_yearmm or die "Failed to create path: $grape_path_yearmm"; }
    # create unique filenames
    my $grape_metadata = substr($stamp,0,6) . "/ISS_" . $stamp . ".tsv";
    my $grape_tree = substr($stamp,0,6) . "/ISS_" . $stamp . ".nwk";
    copy($phantg_metadata,$grape_path . $grape_metadata);
    copy($phantg_tree,$grape_path . $grape_tree);
    # create the html file linking the tree and metadata files
    open my $of, '>', "$phantg_grapetree" or die "Cannot open $phantg_grapetree: $!";
    print $of "<!DOCTYPE html><html><body>\n";
    print $of "<strong><a href=\"https:\/\/irida.iss.it\/grapetree\/?tree=$grape_tree\&metadata=$grape_metadata\">visualizza albero filogenetico<\/a><\/strong>\n";
    print $of "<\/body><\/html>\n";
    close $of;
}