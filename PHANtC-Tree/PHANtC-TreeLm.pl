#!/usr/bin/env perl
## A wrapper script to call reporTree
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;
use DBI;
use Config::Simple;
use JSON::PP;

# Parse arguments
my ($sample_code,
    $metadata_json,
    $phantclm_tree,
    $phantclm_dm,
    $phantclm_grapetree,
    $phantclm_cluster) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my $sample_id = substr($sample_code,2);
my $sample_metadata = readJsonFile();
if ($sample_metadata eq "ND") {
    open my FILEHANDLE, '>', $phantclm_tree and close FILEHANDLE or die "Failed to create file: $!\n";
    open my FILEHANDLE, '>', $phantclm_dm and close FILEHANDLE or die "Failed to create my : $!\n";
    open my FILEHANDLE, '>', $phantclm_grapetree and close FILEHANDLE or die "Failed to create my : $!\n";
    open my FILEHANDLE, '>', $phantclm_cluster and close FILEHANDLE or die "Failed to create my : $!\n";
} else {
    createAllelesFile();
    createMetadataFile();
    runReporTree();
}
exit(0);

# obtain metadata from json file
sub readJsonFile {
    open my $jf, '<', $metadata_json or die "Can't open file $!";
    read $jf, my $json_text, -s $jf;
    chomp($json_text);
    my $metadata = "ND";
    my $json_var = decode_json substr($json_text,1,-1);
    my $contaminated = $json_var->{qc_messages};
    if (index($contaminated, "contaminated") != -1) {
        my $sequence = $json_var->{information_name};
        my $region = $json_var->{region};
        my $country = "Italy";
        my $date_condizioneclinica_origine = getSampleMetadata();
        my $serogroup = $json_var->{serotype_serogroup};
        my $amplicons = $json_var->{serotype_amplicons};
        my $mlst_st = $json_var->{mlst_ST};
        my $mlst_cc = $json_var->{mlst_CC};
        my $mlst_lineage = $json_var->{mlst_lineage};
		$metadata = join("\t", $sequence, $region, $country, $date_condizioneclinica_origine, $serogroup, $amplicons, $mlst_st, $mlst_cc, $mlst_lineage)
	}
    return $metadata;
}

# Obtain metadata from db not present in json file
sub getSampleMetadata {
    my $sql = "select collectionDate,isolate,isolationSource from sample where id=$sample_id";
    my $metadata_db;
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    while (my @row = $sth->fetchrow_array) { 
      $metadata_db = join("\t", @row);
    }
    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
    return $metadata_db;
}

# Obtain allele profiles from db with allele_coverage>79,9% and write to temp file
sub createAllelesFile {
    my $sql = "select allele_strain from mlst_listeria where sample_code='FILE' union
               select allele_strain from mlst_listeria where permille_loci>799";
    # connect to MySQL database
    my %attr2 = ( PrintError=>0, RaiseError=>1);
    my $dbh2 = DBI->connect($dsn,$user,$pwd,\%attr2);
    my $sth2 = $dbh2->prepare($sql);
    $sth2->execute();

    open my $of, '>', "cgMLST.tsv" or die "Cannot open cgMLST.tsv: $!";
    while(my @row = $sth2->fetchrow_array()){
      print $of "$row[0]\n"
    }       
    close $of;

    $sth2->finish();
    # disconnect from the MySQL database
    $dbh2->disconnect();
    return 0;
}

# Obtain metadata from db and write to output file
sub createMetadataFile {
    my $sql = "select Ceppo,Regione,InizioSintomi,CondizioneClinica,Origine,Serogroup,Amplicons,MLST_ST,MLST_CC,MLST_Lineage from v_listeria_opendata";
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "phantclm_metadata.tsv" or die "Cannot open phantclm_metadata.tsv: $!";
    print $if "sequence\tregion\tcountry\tdate\tCondizioneClinica\tOrigine\tSerogroup\tAmplicons\tMLST ST\tMLST CC\tlineage\n";
    print $if "$sample_metadata\n";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[0]\t$row[1]\tItaly\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\n";
    }
    close $if;
    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
}

# Run ReporTree
sub runReporTree {
    # calc distance matrix and minimum spanning tree
    my $result = system("python $scriptdir/reportree.py -a cgMLST.tsv -m phantclm_metadata.tsv --analysis grapetree -thr 4,7,15");
    copy("ReporTree_dist_hamming.tsv", $phantclm_dm);
    copy("ReporTree.nwk", $phantclm_tree);
    copy("ReporTree_partitions.tsv", $phantclm_cluster);
    return 0;
}

sub createGrapeTreeLink {
    # copy phantclm_metadata.tsv & $phantclm_tree to a shared path visible by GrapeTree on the server vizapp.iss.it
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
    copy("phantclm_metadata.tsv",$grape_path . $grape_metadata);
    copy($phantclm_tree,$grape_path . $grape_tree);
    # create the html file linking the tree and metadata files
    my $strUrl = "https://irida.iss.it/grapetree/?tree=$grape_tree&metadata=$grape_metadata";
    open my $of, '>', "$phantclm_grapetree" or die "Cannot open $phantclm_grapetree: $!";
    print $of "<!DOCTYPE html><html><head>";
    print $of "<meta http-equiv=\"refresh\" content=\"0; URL=$strUrl\" />\n";
    print $of "</head><body>\n";
    print $of "<strong><a href=\"$strUrl\">visualizza albero filogenetico<\/a><\/strong>\n";
    print $of "<\/body><\/html>\n";
    close $of;
}