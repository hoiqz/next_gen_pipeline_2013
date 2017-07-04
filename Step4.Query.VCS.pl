#!/usr/bin/perl -w

###########################################################################
#
# We query each variant to see if it is a VCS, and mark as such.
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


if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}

my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;
print "Script: $0\n";
my $TRUE = 0;
my $FALSE = 1;

my $DEBUG = $FALSE; # Set to $FALSE during production

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
my $insubdir = $subdirs{VEP_ANNOTATED};
my $indir = "$subrootdir/$insubdir";
my $inFile = "$indir/$rootname.annotated.by.vep.CAG.vcf"; # At this point, the file is VCF
my $isVCF = $userConfig{VEP_INFILE_FORMAT};


my $inFile_not_CAG = "$indir/$rootname.annotated.by.vep.not.CAG.vcf"; # At this point, the file is VCF


#if ($DEBUG == $TRUE) {
    # I just want to test that my handling of indels are sound, injecting fake insertions/deletions into called variants 
#    print "I am debugging\n";
#    my $debug_addition_cag = "/home/polaris/t1/VIR/injection/INDEL.TEST.CAG.vep";
#    my $debug_addition_not_cag = "/home/polaris/t1/VIR/injection/INDEL.TEST.not.CAG.vep";
#    print "cat $debug_addition_not_cag >> $inFile_not_CAG\n";
#    print "cat $debug_addition_cag >> $inFile\n";
#    `cat $debug_addition_not_cag >> $inFile_not_CAG`;
#    `cat $debug_addition_cag >> $inFile`;
#}



# CREATE OUTPUT 
my $outsubdir = $subdirs{MARK_VAR_CALLERS};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) {
    system("mkdir -p $outdir");
}
my $outfile_CAG = "$outdir/$rootname.marked.CAG.vcs";

my $outfile_not_CAG = "$outdir/$rootname.marked.not.CAG.vcs";


my $journal_file = "$outdir/$rootname.journal.references.CAG.txt";
my $generic_genes_file = "$outdir/$rootname.generic.genes.CAG.txt";

my $journal_file_not_CAG = "$outdir/$rootname.journal.references.not.CAG.txt";
my $generic_genes_file_not_CAG = "$outdir/$rootname.generic.genes.not.CAG.txt";


# Database
my $database = $references{REFLEX_VCS_DATABASE};
my $table = $references{REFLEX_VCS_DATABASE_TABLE};
my $gene_database = $references{CLINICALLY_ACTIONABLE_GENES_DB};
my $gene_table = $references{CLINICALLY_ACTIONABLE_GENES_TABLE};    


my ($generic_genes_href, $journal_references_href) = &markIfVCS($database, $table, $disease, $inFile, "vcf", $outfile_CAG, $gene_database, $gene_table);
&writeFile($generic_genes_href, $generic_genes_file);
&writeFile($journal_references_href, $journal_file);



print "Completed: $outfile_CAG\n";



my ($generic_genes_not_CAG_href, $journal_references_not_CAG_href) = 
    &markIfVCS($database, $table, $disease, $inFile_not_CAG, "vcf", $outfile_not_CAG, $gene_database, $gene_table);
&writeFile($generic_genes_not_CAG_href, $generic_genes_file_not_CAG);
&writeFile($journal_references_not_CAG_href, $journal_file_not_CAG);

print "Completed: $outfile_not_CAG\n";


sub writeFile() {
    my ($href, $outfile) = @_;
    print "Writing out to $outfile\n";
    my %hash = %{$href};
    my @data = keys(%hash);
    open(OUT, ">$outfile") || die "Unable to open $outfile for writing.\n";
    foreach my $row (@data) {
	chomp $row;
	print OUT "$row\n";
    }#end foreach
    close(OUT);
}#end writeFile



__END__
