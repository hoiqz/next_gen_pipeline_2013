#!/usr/bin/perl -w

###########################################################################
#
# Prepare for Validation Statistics
# 
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.10.25
#
###########################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
use DBI;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.vcs.checker.pl";
require "$pwd/polaris.independent.validation.checker.pl";

if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}

my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;
print "Script: $0\n";
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
my $insubdir = $subdirs{VARSCAN_DIR};
my $indir = "$subrootdir/$insubdir";
my $infile = "$indir/$rootname.all.variants";

my $isVCF = $userConfig{VEP_INFILE_FORMAT};

# OUTFILE
my $outsubdir = $subdirs{VALIDATION_STATS_DIR};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $intermediate_file = "$outdir/$rootname.validation.mark.vcs";
my $final_file = "$outdir/$rootname.for.validation";

# Now, identify which are VCS
my $database = $references{REFLEX_VCS_DATABASE};
my $table = $references{REFLEX_VCS_DATABASE_TABLE};

my $gene_database = $references{CLINICALLY_ACTIONABLE_GENES_DB};
my $gene_table = $references{CLINICALLY_ACTIONABLE_GENES_TABLE};

my $iv_database = $references{REFLEX_INDEPENDENTLY_VALIDATED_DB};
my $iv_table = $references{REFLEX_INDEPENDENTLY_VALIDATED_TABLE};

my ($generic_genes_href, $journal_references_href) = &markIfVCS($database, $table, $disease, $infile, $isVCF, $intermediate_file, $gene_database, $gene_table);

# Annotate independent validation
&annotate_independent_validation($iv_database, $iv_table, $disease, $intermediate_file, $final_file);

print "Completed: $final_file\n";






__END__

