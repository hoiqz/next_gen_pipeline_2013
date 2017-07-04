#!/usr/bin/perl -w

################################################################################
#
# For each variant, we check if they have been independently validated or not
# 
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.10.28
#
################################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
use DBI;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.independent.validation.checker.pl";

if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}
print "Script: $0\n";
my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;

my $TRUE = 0;
my $FALSE = 1;

my $INDP_VALID = "IS_INDP_VALIDATED";
my $NO_INDP_VALID = "NOT_INDP_VALIDATED";

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
my $insubdir = $subdirs{SELECTED_VEP_ANNOTATION};
my $indir = "$subrootdir/$insubdir";
my $infile = "$indir/$rootname.CAG.selected.annotations";


my $infile_not_CAG = "$indir/$rootname.not.CAG.selected.annotations";


# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{CHECK_INDEPENDENT_VALIDATION_DIR};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $outfile = "$outdir/$rootname.indp.validate.CAG.annotated";

my $outfile_not_CAG = "$outdir/$rootname.indp.validate.not.CAG.annotated";


# DATABASE INFORMATION
my $database = $references{REFLEX_INDEPENDENTLY_VALIDATED_DB};
my $table = $references{REFLEX_INDEPENDENTLY_VALIDATED_TABLE};

&annotate_independent_validation($database, $table, $disease, $infile, $outfile);

&annotate_independent_validation($database, $table, $disease, $infile_not_CAG, $outfile_not_CAG);

print "Completed: $outfile\n";
print "Completed: $outfile_not_CAG\n";

__END__

