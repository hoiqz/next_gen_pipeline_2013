#!/usr/bin/perl -w

###########################################################################
#
# We call indels using dindel (http://www.sanger.ac.uk/resources/software/dindel/)
#
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.12.18
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

#if (scalar(@ARGV) != 3) {
#    die "Usage: perl $0 <bamfile> <name> <outdir>\n";
#}

#my ($bamfile, $name, $outdir) = @ARGV;
my $subdir = $subdirs{PILEUP_DIR};
my $dir = "$subrootdir/$subdir";
my $bamfile = "$dir/$rootname" . "_reAligned_reCal.bam";

my $name = $rootname;

# Create output directory
my $outsubdir = $subdirs{INDEL_DIR};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $intermediate = "dindel-intermediate";
my $intermediate_dir = "$outdir/$intermediate";
system("mkdir -p $intermediate_dir");
my $outfile = "$outdir/$name" . ".vcf";

# Get executables
my $python_exec = $softwareTools{PYTHON_EXEC};
my $dindel = $softwareTools{DINDEL_EXEC};
my $makeWindows = $softwareTools{DINDEL_MAKEWIN};
my $mergeScript = $softwareTools{DINDEL_MERGE};
my $genome_dir = $references{GENOME_REF_DIR};
my $build = $userConfig{REF_BUILD};
my $ref = $genome_dir . "/" . $build . "/" . $build . ".fa";



# Step 1
my $step1_cmd = "$dindel --analysis getCIGARindels --bamFile $bamfile --outputFile $intermediate_dir/$name --ref $ref";
print "$step1_cmd\n";
my $result = system($step1_cmd);
print "Step 1: $intermediate_dir/$name.variants and $intermediate_dir/$name.libraries.txt\n";

#Step 2: $makeWindows --inputVarFile $varfile --windowFilePrefix $name.realign_windows --numWindowsPerFile 1000
my $varfile = "$intermediate_dir/$name.variants.txt";
my $prefix = "$intermediate_dir/$name.realign_windows";
my $step2_cmd = "$python_exec $makeWindows --inputVarFile $varfile --windowFilePrefix $prefix --numWindowsPerFile 1000";
print "$step2_cmd\n";
$result = system($step2_cmd);

#Step 3: Get number of realign_windows files
my $num = 0;
opendir(DIR, $intermediate_dir) || die "Unable to open $intermediate_dir\n";
while(my $file = readdir(DIR)) {
    if ($file =~ /$name\.realign_windows\.(\d+)\.txt/) {
	if ($num < $1) { $num = $1; }
    }
}#end while
closedir(DIR);
print "Number of re-aligned windows: $num\n";


#Step 3: dindel --analysis indels --doDiploid --bamFile sample.bam --ref ref.fa \
#    --varFile sample.realign_windows.2.txt \
#    --libFile sample.dindel_output.libraries.txt \
#    --outputFile sample.dindel_stage2_output_windows.2
my @childs = ();
my $prefix_name = "$intermediate_dir/$name";
for ( my $count = 1; $count <= $num; $count++) {
    my $pid = fork();
    if ($pid != 0) {
        # print "pid is $pid, parent $$\n";
        push(@childs, $pid);         # parent
    } elsif ($pid == 0) {
	&runStep3($count, $prefix_name); # child
	exit 0;
    } else {
	die "Could not fork: $!\n";
    }
}
 
foreach (@childs) {
    my $tmp = waitpid($_, 0);
    print "Done with pid $tmp\n"; 
}


#Step 4: mergeOutputDiploid.py --inputFiles sample.dindel_stage2_outputfiles.txt --outputFile variantCalls.VCF --ref ref.fa
# Now, we need to create a list containing *glf.txt files 

my @glf = ();
opendir(DIR, $intermediate_dir) || die "Unable to open $intermediate_dir\n";
while(my $file = readdir(DIR)) {
    if ($file =~ /\.glf\.txt/) { 
	push @glf, "$intermediate_dir/$file"; 
    }
}#end while
closedir(DIR);
my $glf_list = "$intermediate_dir/glf.txt";
open(FILE, ">$glf_list") || die "Unable to open $glf_list\n";
foreach my $glf_file (@glf) {
    chomp $glf_file;
    print FILE "$glf_file\n";
}
close(FILE);

my $step4_cmd = "$python_exec $mergeScript --inputFiles $glf_list --outputFile $outfile --ref $ref";
print "$step4_cmd\n";
$result = system($step4_cmd);

print "Completed: $outfile\n";

my $filtered_outfile = "$outfile" . ".filtered"; # VEP has trouble with chromosomes like chr17_gl000204_random
my @chrs = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrx chry chrM chrm);
my %chromosomes = ();
foreach my $c (@chrs) {
    $chromosomes{$c} = 0;
}
open(OUT, ">$filtered_outfile") || die "Unable to open for writing: $filtered_outfile\n";
open(IN, "<$outfile") || die "Unable to open $outfile for reading.\n";
while(my $row = <IN>) {
    chomp $row;
    if ($row =~ /^\#/) { print OUT "$row\n"; }
    else {
	my ($chr, @others) = split(/\t/, $row);
	if (defined($chromosomes{$chr})) {
	    print OUT "$row\n";
	}
    }
}
close(IN);
close(OUT);

# Now we annotate with VEP
my $perl = $softwareTools{PERL_EXEC};
my $vep_script = $softwareTools{VEP_PERLSCRIPT};
my $vep_lib = $softwareTools{VEP_LIB};

my $annotated_outfile = "$outdir/$rootname.indels.from.dindel.annotated.vcf";


my $vep_cmd = "$perl $vep_script --cache -i " . $filtered_outfile . " --format vcf " . # The dindel output file is a VCF file
    "-coding_only --no_intergenic " . # This focuses the annotations in coding regions only, and makes the file much smaller
    "-o " . $annotated_outfile . " --vcf " .  # Output
    "--offline -dir $vep_lib --force_overwrite " . # I/O 
    "--species homo_sapiens " . # Build information
    "--ccds --xref_refseq --gmaf --maf_1kg  --numbers --polyphen b"; # Options

print "$vep_cmd\n";
$result = system($vep_cmd);
print "Completed: $annotated_outfile\n";




sub runStep3() {
    my ($i, $prefix) = @_;
    my $cmd = "$dindel --analysis indels --doDiploid --bamFile $bamfile --ref $ref " . 
	"--varFile $prefix.realign_windows.$i.txt " .
	"--libFile $prefix.libraries.txt " . 
	"--outputFile $prefix.dindel_stage2_output_windows.$i";
    print "$cmd\n";
    $result += system($cmd);
}#end runStep3



__END__
