#!/Usr/bin/perl -w

################################################################################
#
# We take all variants found and annotate using ANNOVAR,
# then parse the information to get what we want.
# 
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.11.04
#
################################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
#use Class::Struct;
use DBI;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.objects.pl";
require "$pwd/polaris.gather.data.for.recording.pl";



if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}
print "Script: $0\n";
my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;

my $TRUE = 0;
my $FALSE = 1;

my ($sw_href, $subdir_href, $user_href, $references_href) = &getConfigurations($swconfigfile, $subdirconfigfile, $userconfigfile, $referencesfile);
my %softwareTools = %{$sw_href};
my %subdirs = %{$subdir_href};
my %userConfig = %{$user_href};
my %references = %{$references_href};

# Get root directory
my ($rootdir, $hospital, $patientid, $disease) = &getRootDirInfo($user_href);

my $rootname = "$hospital.$patientid.$disease";
my $subrootdir = "$rootdir/$rootname";

# INPUT INFORMATION
my $vep_insubdir = $subdirs{SELECTED_VEP_ANNOTATION};
my $vep_indir = "$subrootdir/$vep_insubdir";

my $vep_CAG_infile = "$vep_indir/$rootname.CAG.selected.annotations"; 
my $vep_not_CAG_infile = "$vep_indir/$rootname.not.CAG.selected.annotations"; 
my $vep_infile = "$vep_indir/$rootname.all.selected.annotations"; 
&concatenateVEP($vep_CAG_infile, $vep_not_CAG_infile, $vep_infile);


my $varscan_insubdir = $subdirs{VARSCAN_DIR};
my $varscan_indir = "$subrootdir/$varscan_insubdir";
my $varscan_CAG_infile = "$varscan_indir/$rootname.CAG.variants";

my $varscan_not_CAG_infile = "$varscan_indir/$rootname.not.CAG.variants";
my $varscan_infile = "$varscan_indir/$rootname.all.variants";


my $outsubdir = $subdirs{DEPOSIT_FOR_RECORDING};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir $outdir"); }
my $outputfile = "$rootname.for.recording";
my $outfile = "$outdir/$outputfile";
&concatenateVarScan($varscan_CAG_infile, $varscan_not_CAG_infile, $varscan_infile);





&gatherData($outdir, $rootname, $vep_infile, $varscan_infile, $disease, $hospital, $patientid, $outfile);



my $patientConsent = "YES";
$patientConsent = $userConfig{CONSENT_TO_STORE};
if ($patientConsent && $patientConsent =~ /NO/i) {
    # Then the patient has NOT consented to have his results stored in Detected Variant database
    print "Completed: Patient has elected NOT to have his sequenced results stored in Detected Variants Database\n";
    exit(0);
}


# CREATE OUTPUT DIRECTORY
my $ROOTDIR = $references{DETECTED_VARIANTS_CRON_DIR};
if (! -d $ROOTDIR) { system("mkdir -p $ROOTDIR"); }

my $output_sub_dir = &createOutputSubDirectory($ROOTDIR, $hospital, $patientid, $disease);

# Make a sub directory for configuration files
my $config_subdir = "$output_sub_dir/config";
system("mkdir -p $config_subdir");
system("cp $swconfigfile $config_subdir");
system("cp $subdirconfigfile $config_subdir");
system("cp $referencesfile $config_subdir");
system("cp $userconfigfile $config_subdir");

# We now need the file name of the configuration files
my $config_sw = &getFileName($swconfigfile);
my $config_subdirectories = &getFileName($subdirconfigfile);
my $config_ref = &getFileName($referencesfile);
my $config_user = &getFileName($userconfigfile);

# Now we create a meta file
my $metafile = "$config_subdir/META.INFO";
open(META, ">$metafile") || die "Unable to open $metafile for writing.\n";
print META "POLARIS_SOFTWARE:$config_sw\n";
print META "POLARIS_SUBDIRS:$config_subdirectories\n";
print META "POLARIS_REFERENCE:$config_ref\n";
print META "POLARIS_USERCONF:$config_user\n";

my $time = time();
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
$year += 1900; $mon += 1;
print META "DATE_OF_RUN:$year,$mon,$mday,$hour,$min,$sec,$time\n";

close(META);

print "$config_subdir/META.INFO\n";

my $ofile = $output_sub_dir . "/" . $references{DETECTED_VARIANTS_CRON_INFO_FILE};

print "Copying file $outfile\n";
print "to $ofile\n";
my $result = system("cp $outfile $ofile");
if ($result == 0) {
    print "Completed: $ofile\n";
} else {
    print "FAILED cp $outfile $ofile\n";
}



### SUB-ROUTINE ########################################################

sub createOutputSubDirectory() {
    my ($dir, $id1, $id2, $id3) = @_;
    my $time = time();
    my $directory = "$dir/$id1.$id2.$id3.$time";
    system("mkdir -p $directory");
    print "$directory\n";
    return $directory;
}#end createOutputSubDirectory

sub getFileName() {
    my $file = $_[0];
    print "Original: $file\n";
    my @suffix = ();
    my $basename = basename($file,@suffix);
    print "New: $basename\n";
    return $basename;
}#end getFileName

sub concatenateVEP() {
    my ($vep_CAG_infile, $vep_not_CAG_infile, $vep_infile) = @_;
    system("cp $vep_CAG_infile $vep_infile");
    # This is correct, it is an infile, but we are appending to it.
    open(OUT, ">>$vep_infile") || die "Unable to open $vep_infile for writing.\n";
    open(IN, "<$vep_not_CAG_infile") || die "Unable to open $vep_not_CAG_infile for reading.\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row !~ /\#/) {
	    print OUT "$row\n";
	}
    }#end while
    close(IN);
    close(OUT);

}#end concatenateVEP

sub concatenateVarScan() {
    my ($varscan_CAG_infile, $varscan_not_CAG_infile, $varscan_infile) = @_;
    system("cp $varscan_CAG_infile $varscan_infile");
    open(OUT, ">>$varscan_infile") || die "Unable to open $varscan_infile\n";
    open(IN, "<$varscan_not_CAG_infile") || die "Unable to open $varscan_not_CAG_infile\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row !~ /^Chrom/) {
	    print OUT "$row\n";
	}
    }
    close(IN);
    close(OUT);
}#end concatenateVarScan


__END__

