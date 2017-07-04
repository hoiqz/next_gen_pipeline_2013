#!/usr/bin/perl -w

###########################################################################
#
# We use VEP to annotate each variant found in Clinically Actionable Genes
# but that are NOT.VCS
#
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.10.28
#
###########################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
use DBI;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";

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
my $infpath = "$indir/$rootname.CAG.variants";

my $infpath_not_CAG = "$indir/$rootname.not.CAG.variants";



my $infile_format = lc($userConfig{VEP_INFILE_FORMAT}); # pileup or vcf, must be in lowercase for VEP

# CREATE OUTPUT FILE
my $outsubdir = $subdirs{VEP_ANNOTATED};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }

my $outfile = "$rootname.annotated.by.vep.CAG.vcf";
my $outfpath = "$outdir/$outfile";

my $outfile_not_CAG = "$rootname.annotated.by.vep.not.CAG.vcf";
my $outfpath_not_CAG = "$outdir/$outfile_not_CAG";

my $perl = $softwareTools{PERL_EXEC};
my $vep_script = $softwareTools{VEP_PERLSCRIPT};
my $vep_lib = $softwareTools{VEP_LIB};

# --force_overwrite = if output file exists, we force script to overwrite.
# --ccds = adds CCDS
# --xref_refseq = this gives us NM_XXXX so that we can identify the variant IF Disease Champions gives us NM_ instead
# --numbers = Add exon/intron number affected by mutation
# --gmaf = global minor allele freq from 1000 Genomes Phase 1
# --maf_1kg = Add MAF from continental populations (AFR,AMR,ASN,EUR) of 1000 Genomes Phase 1
# --no_intergenic = filters out annotations that are intergenic
# -coding_only = includes only coding annotations

# OTHER FILTERING THAT WE MIGHT CONSIDER IN THE FUTURE
# --filter_common


my $cmd = "$perl $vep_script --cache -i " . $infpath . " --format $infile_format " . # Input
    "-coding_only --no_intergenic " . # This focuses the annotations in coding regions only, and makes the file much smaller
    "-o " . $outfpath . " --vcf " .  # Output
    "--offline -dir $vep_lib --force_overwrite " . # I/O 
    "--species homo_sapiens " . # Build information
    "--ccds --xref_refseq --gmaf --maf_1kg  --numbers --polyphen b"; # Options


print "$cmd\n";
my $results = system($cmd);
if ($results == 0) {
    print "Completed: $outfpath\n";
} else {
    print "Failed to annotate with VEP: $outfpath\n";
    print "$cmd\n";
}


my $cmd_not_CAG = "$perl $vep_script --cache -i " . $infpath_not_CAG . " --format $infile_format " . # Input
    "-coding_only --no_intergenic " . # This focuses the annotations in coding regions only, and makes the file much smaller
    "-o " . $outfpath_not_CAG . " --vcf " .  # Output
    "--offline -dir $vep_lib --force_overwrite " . # I/O 
    "--species homo_sapiens " . # Build information
    "--ccds --xref_refseq --gmaf --maf_1kg  --numbers --polyphen b"; # Options


print "$cmd_not_CAG\n";
my $results_not_CAG = system($cmd_not_CAG);
if ($results_not_CAG == 0) {
    print "Completed: $outfpath_not_CAG\n";
} else {
    print "Failed to annotate with VEP: $outfpath_not_CAG\n";
    print "$cmd_not_CAG\n";
}



__END__

=pod 
# From Variant Effect Predictor website, I cannot find any way to "force" it to use build 37, but I'm not adding this because of the warning below:
--db_version [number]
--db Force the script to connect to a specific version of the Ensembl databases. Not recommended as there will usually be conflicts between software and database versions. Not used by default
-fork is unstable and fails sometimes, I just removed it. It will be slower but at least it won't fail.
http://asia.ensembl.org/info/docs/variation/vep/vep_script.html (Forking part)
=cut

