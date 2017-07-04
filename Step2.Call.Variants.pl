#!/usr/bin/perl -w

###########################################################################
#
# We now call variants in clinically actionable genes for disease
# and query VCS database
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
my $indir = $subdirs{LIMITED_DIR};
my $inFile = "$subrootdir/$indir/$rootname.limited.to.clinically.actionable.genes";

my $inFile_not_CAG = "$subrootdir/$indir/$rootname.not.in.clinically.actionable.genes";


# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{VARSCAN_DIR};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) {
    system("mkdir -p $outdir");
}
my $outSNPfile = "$outdir/$rootname.snps.from.variants";
my $outINDELfile = "$outdir/$rootname.indels.from.varscan.variants";
my $out_pre_filtering_file = "$outdir/$rootname.pre.filtering.variants";
my $outfile = "$outdir/$rootname.CAG.variants";


my $outSNPfile_not_CAG = "$outdir/$rootname.snps.not.CAG.from.variants";
my $outINDELfile_not_CAG = "$outdir/$rootname.indels.from.varscan.not.CAG.variants";
my $out_pre_filtering_file_not_CAG = "$outdir/$rootname.pre.filtering.not.CAG.variants";
my $outfile_not_CAG = "$outdir/$rootname.not.CAG.variants";

my $outfile_all_variants = "$outdir/$rootname.reflex.varscan";

if ($userConfig{VEP_INFILE_FORMAT} =~ /vcf/i) {
    system("cp $inFile $outfile"); # The vcf file limited to clinically actionable genes is the output file in this case, since we're not calling variants using VarScan2
} else {
    # Get varscan tools, etc.
    my $java = $softwareTools{JAVA_EXEC};
    my $varscanJar = $softwareTools{VARSCAN_JAR};
    my $indelOpts = $softwareTools{VARSCAN_INDEL_OPTS};
    my $snpOpts = $softwareTools{VARSCAN_SNP_OPTS};

    my $min_forward_reads = $softwareTools{MIN_ALT_READS_FORWARD_STRAND};
    my $min_reverse_reads = $softwareTools{MIN_ALT_READS_REVERSE_STRAND};

    # Clinically Actionable Genes
    &varscanIt($inFile, $outSNPfile, $java, $varscanJar, "pileup2snp", $snpOpts);
    &varscanIt($inFile, $outINDELfile, $java, $varscanJar, "pileup2indel", $indelOpts);
    &joinFiles($outSNPfile, $outINDELfile, $out_pre_filtering_file);
    &filterByStrandSupport($out_pre_filtering_file, $outfile, $min_forward_reads, $min_reverse_reads);


    # NOT in Clinically Actionable Genes (for Discovery Report)
    &varscanIt($inFile_not_CAG, $outSNPfile_not_CAG, $java, $varscanJar, "pileup2snp", $snpOpts);
    &varscanIt($inFile_not_CAG, $outINDELfile_not_CAG, $java, $varscanJar, "pileup2indel", $indelOpts);
    &joinFiles($outSNPfile_not_CAG, $outINDELfile_not_CAG, $out_pre_filtering_file_not_CAG);
    &filterByStrandSupport($out_pre_filtering_file_not_CAG, $outfile_not_CAG, $min_forward_reads, $min_reverse_reads);


}


if ($DEBUG == $TRUE) {
    `cat /home/polaris/t1/VIR/injection/INDEL.TEST.CAG.varscan >> $outfile`;
    `cat /home/polaris/t1/VIR/injection/INDEL.TEST.not.CAG.varscan >> $outfile_not_CAG`;
}#end if ($DEBUG == $TRUE)


`cat $outfile $outfile_not_CAG > $outfile_all_variants`; # We concatenate both CAG and non-CAG for Step7.Generate.Clinical.Report.pl


print "Completed: $outfile\n";
print "Completed: $outfile_not_CAG\n";


## SUB-ROUTINES


sub filterByStrandSupport() {

=pod
From Simeen's Methods section:
The variant information provided in this report was generated using sequenced files generated from MiSeq and processed with reference to hg19 as reference genome. Variants are called if they have an minimum read coverage of 100, variant frequency of at least 0.05, alternate forward strand detected with at least 10 reads, alternative reverse strand detected with at least 10 reads, and minimum base-call quality of 30.

The minimum read coverage (100), variant frequency (0.05), and minimum base-call quality (30) were handled by VARSCAN_SNP_OPTS.
This sub-routine handles the alternate forward strand (10 reads) and alternative reverse strand (10 reads)
=cut

    my ($infile, $outfile, $min_alt_forward_reads, $min_alt_reverse_reads) = @_;
    # These are VarScan files
    open(OUT, ">$outfile") || die "Unable to open $outfile for writing.\n";
    open(IN, "<$infile") || die "Unable to open $infile for reading.\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row =~ /^Chrom/) {
	    print OUT "$row\n";
	} else {
	    my ($Chrom, $Position, $Ref, $Cons, $Reads1, $Reads2, $VarFreq, 
		$Strands1, $Strands2, $Qual1, $Qual2, $Pvalue, $MapQual1, $MapQual2,
		$Reads1Plus, $Reads1Minus, $Reads2Plus, $Reads2Minus, $VarAllele) = split(/\t/, $row);
	    
	    if (($min_alt_forward_reads <= $Reads2Plus) && ($min_alt_reverse_reads <= $Reads2Minus)) {
		print OUT "$row\n";
	    }
	}
    }#end while
    close(IN);
    close(OUT);

}#end filterByStrandSupport

sub joinFiles() {
    my ($snp, $indel, $outfile) = @_;
    system("cp $snp $outfile");
    `grep -v "Chrom" $indel >> $outfile`;
}#end joinFiles

sub varscanIt() {
    my ($infile, $outfile, $java, $varscan_jar, $tool, $opts) = @_;

    my $cmd = "$java -jar $varscan_jar $tool $infile $opts > $outfile";
    print "$cmd\n";
    my $results = system($cmd);
}#end handleSNPs




__END__

