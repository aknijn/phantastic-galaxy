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
    $phantcec_tree,
    $phantcec_dm,
    $phantcec_grapetree,
    $phantcec_cluster) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
my $cfg = new Config::Simple("$scriptdir/../phantastic.conf");
my $dsn = $cfg->param('db.dsn');
my $user = $cfg->param('db.user');
my $pwd = $cfg->param('db.password');
my (undef, $sample_id) = split('_', $sample_code);
my ($antigen_o, $sample_metadata) = readJsonFile();
if ($sample_metadata eq "ND") {
    open FILEHANDLE, '>', $phantcec_tree and close FILEHANDLE or die "Failed to create file: $!\n";
    open FILEHANDLE, '>', $phantcec_dm and close FILEHANDLE or die "Failed to create file: $!\n";
    open FILEHANDLE, '>', $phantcec_grapetree and close FILEHANDLE or die "Failed to create file: $!\n";
    open FILEHANDLE, '>', $phantcec_cluster and close FILEHANDLE or die "Failed to create file: $!\n";
} else {
    createAllelesFile();
    createMetadataFile();
    runReporTree();
    createGrapeTreeLink();
}
exit(0);

sub getLatLong {
	# Italian coordinates
	my $inRegion = shift;
	my $lat = "42.833333";
    my $long = "12.833333";
    if ($inRegion eq "Piemonte") { $lat = "45.066667"; $long = "7.700000"; }
    if ($inRegion eq "Valle d'Aosta") { $lat = "45.737222"; $long = "7.320556"; }
    if ($inRegion eq "Lombardia") { $lat = "45.464161"; $long = "9.190336"; }
    if ($inRegion eq "Trentino-Alto Adige") { $lat = "11.116667"; $long = "46.066667"; }
    if ($inRegion eq "Veneto") { $lat = "45.439722"; $long = "12.331944"; }
    if ($inRegion eq "Friuli-Venezia Giulia") { $lat = "13.804167"; $long = "45.636111"; }
    if ($inRegion eq "Liguria") { $lat = "44.411156"; $long = "8.932661"; }
    if ($inRegion eq "Emilia-Romagna") { $lat = "44.493889"; $long = "11.342778"; }
    if ($inRegion eq "Toscana") { $lat = "43.771389"; $long = "11.254167"; }
    if ($inRegion eq "Umbria") { $lat = "43.112100"; $long = "12.388800"; }
    if ($inRegion eq "Marche") { $lat = "43.616667"; $long = "13.516667"; }
    if ($inRegion eq "Lazio") { $lat = "41.893056"; $long = "12.482778"; }
    if ($inRegion eq "Abruzzo") { $lat = "42.354008"; $long = "13.391992"; }
    if ($inRegion eq "Molise") { $lat = "41.561000"; $long = "14.668400"; }
    if ($inRegion eq "Campania") { $lat = "40.833333"; $long = "14.250000"; }
    if ($inRegion eq "Puglia") { $lat = "41.125278"; $long = "16.866667"; }
    if ($inRegion eq "Basilicata") { $lat = "40.633333"; $long = "15.800000"; }
    if ($inRegion eq "Calabria") { $lat = "38.910000"; $long = "16.587500"; }
    if ($inRegion eq "Sicilia") { $lat = "38.115556"; $long = "13.361389"; }
    if ($inRegion eq "Sardegna") { $lat = "39.216667"; $long = "9.116667"; }
	return ($lat, $long);
}

# obtain metadata from json file
sub readJsonFile {
    open my $jf, '<', $metadata_json or die "Can't open file $!";
    read $jf, my $json_text, -s $jf;
    chomp($json_text);
    my $metadata = "ND";
    my $antigen_o = "ND";
    my $json_var = decode_json substr($json_text,1,-1);
    my $contaminated = $json_var->{qc_messages};
    if (index($contaminated, "contaminated") == -1) {
        my $sequence = $json_var->{information_name};
        my $region = $json_var->{region};
        my $country = "Italy";
        my $date_campione_condizioneclinica_origine = getSampleMetadata();
        $antigen_o = $json_var->{serotype_o};
        my $antigen_h = $json_var->{serotype_h};
        my $mlst_st = $json_var->{mlst_ST};
        my $stx1 = $json_var->{virulotype_stx1};
        my $stx2 = $json_var->{virulotype_stx2};
        my $stxsubtype = $json_var->{shigatoxin_subtype};
        my $eae = $json_var->{virulotype_eae};
        my $ehxa = $json_var->{virulotype_ehxa};
		my ($latitude, $longitude) = getLatLong($region);
        $metadata = join("\t", $sequence, $region, $country, $date_campione_condizioneclinica_origine, $antigen_o, $antigen_h, $mlst_st, $stx1, $stx2, $stxsubtype, $eae, $ehxa, $latitude, $longitude)
    }
    return ($antigen_o, $metadata);
}

# Obtain metadata from db not present in json file
sub getSampleMetadata {
    my $sql = "select collectionDate,sampleName,IFNULL(isolate,''),IFNULL(isolationSource,'') from sample where id=$sample_id";
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
    my $sql;
    # If the serogroup is known, retrieve allele profiles only from those samples, if O? then from all
    if ($antigen_o eq 'O?') {
        $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
               select allele_strain from mlst_ecoli where permille_loci>799";
    }
    else
    {
        $sql = "select allele_strain from mlst_ecoli where sample_code='FILE' union
               select allele_strain from mlst_ecoli where sample_code='$sample_code' union
               select allele_strain from v_mlst_ecoli where permille_loci>799 and Antigen_O = '$antigen_o'";
    }
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
    my $sql = "select Ceppo,Regione,InizioSintomi,CondizioneClinica,Origine,Sequence,Antigen_O,Antigen_H,MLST,stx1,stx2,stxsub,eae,ehxa,latitude,longitude from v_ecoli_opendata where Antigen_O ='$antigen_o'";
    # connect to MySQL database
    my %attr = ( PrintError=>0, RaiseError=>1);
    my $dbh = DBI->connect($dsn,$user,$pwd,\%attr);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    open my $if, '>', "phantcec_metadata.tsv" or die "Cannot open phantcec_metadata.tsv: $!";
    print $if "sequence\tregion\tcountry\tdate\tCampione\tCondizioneClinica\tOrigine\tAntigen O\tAntigen H\tMLST ST\tstx1\tstx2\tstx subtype\teae\tehxA\tlatitude\tlongitude\n";
    print $if "$sample_metadata\n";
    while (my @row = $sth->fetchrow_array) { 
      print $if "$row[5]\t$row[1]\tItaly\t$row[2]\t$row[0]\t$row[3]\t$row[4]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\t$row[14]\t$row[15]\n";
    }
    close $if;
    $sth->finish();
    # disconnect from the MySQL database
    $dbh->disconnect();
}

# Run ReporTree
sub runReporTree {
    # calc distance matrix and minimum spanning tree
    my $result = system("python $scriptdir/reportree.py -a cgMLST.tsv -m phantcec_metadata.tsv --analysis grapetree -thr 4,7,15");
    copy("ReporTree_dist_hamming.tsv", $phantcec_dm);
    copy("ReporTree.nwk", $phantcec_tree);
    copy("ReporTree_partitions.tsv", $phantcec_cluster);
    return 0;
}

sub createGrapeTreeLink {
    # copy phantcec_metadata.tsv & $phantcec_tree to a shared path visible by GrapeTree on the server vizapp.iss.it
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
    copy("phantcec_metadata.tsv",$grape_path . $grape_metadata);
    copy($phantcec_tree,$grape_path . $grape_tree);
    # create the html file linking the tree and metadata files
    my $strUrl = "https://irida.iss.it/spread/?tree=spread/$grape_tree&metadata=spread/$grape_metadata";
    open my $of, '>', "$phantcec_grapetree" or die "Cannot open $phantcec_grapetree: $!";
    print $of "<!DOCTYPE html><html><head>";
    print $of "<meta http-equiv=\"refresh\" content=\"0; URL=$strUrl\" />\n";
    print $of "</head><body>\n";
    print $of "<strong><a href=\"$strUrl\">visualizza albero filogenetico<\/a><\/strong>\n";
    print $of "<\/body><\/html>\n";
    close $of;
}