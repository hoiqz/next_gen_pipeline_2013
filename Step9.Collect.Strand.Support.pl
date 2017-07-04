#!/usr/bin/perl -w

################################################################################
#
# For each VCS not found in the sample, we compute the strand support
# from the pileup file.
# 
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.11.05
#
################################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
use Class::Struct;
use DBI;


my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.objects.pl";
require "$pwd/polaris.strand.support.pl";

if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}
print "Script: $0\n";
my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;

my $TRUE = 0;
my $FALSE = 1;
my $DEBUG = $FALSE; # Set to false during production
my $VCS = "VCS"; # Because during implementation, we do not have VCS variants found, so we're using NOT.VCS, this must be changed to VCS during production.


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
# We want to get the pileUp file
my $insubdir = $subdirs{PILEUP_DIR};
my $indir = "$subrootdir/$insubdir";
#my $pileup_file = "$rootname" . "_filtered.txt";
my $pileup_file = "$rootname" . "_pileUp";
my $pileup_fpath = "$indir/$pileup_file";

my $found_subdir = $subdirs{SELECTED_VEP_ANNOTATION};
my $found_dir = "$subrootdir/$found_subdir";
my $foundVariantsFile = "$found_dir/$rootname.CAG.selected.annotations";

# VCS Database
my $database = $references{POLARIS_VALIDATED_DB};
my $table = $references{POLARIS_VALIDATED_TABLE};


# CREATE OUTPUT DIRECTORY
my $variant_not_found_subdir = $subdirs{VCS_NOT_FOUND_DIR};
my $outdir = "$subrootdir/$variant_not_found_subdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $not_found_outfile = "$outdir/$rootname.vcs.not.found.coverage";
my $found_outfile = "$outdir/$rootname.vcs.found.coverage";

# Step 1: Get list of VCS found
my $vcs_found_href = &getVCSFound($foundVariantsFile);


# TEST
#my %test = %{$vcs_found_href};
#my @keys = keys(%test);
#foreach my $e (@keys) { print "$e\n"; }
# END TEST


# Step 2: Get list of VCS not found using VCS database
my ($vcs_found_aref, $vcs_not_found_aref) = &splitVCSFoundAndNotFound($vcs_found_href, $disease, $database, $table);


my @notfoundtest = @{$vcs_not_found_aref};

# TEST
#foreach my $e (@notfoundtest) {
#    my $line = $e->chr . "\t" . $e->ref . "\t" . $e->pos . "\t" . $e->alt;
#    print "NOT FOUND: $line\n";
#}
# END TEST
my @foundtest = @{$vcs_found_aref};
if ($DEBUG == $TRUE) {
    my $v1 = VARIANT_STRAND_SUPPORT->new();
    $v1->chr("chr2");
    $v1->pos(29416572);
    $v1->ref("T");
    $v1->alt("C");
    push @foundtest, $v1; # I just want to test it.

    my $v2 = VARIANT_STRAND_SUPPORT->new();
    $v2->chr("chr2");
    $v2->pos(29444095);
    $v2->ref("C");
    $v2->alt("T");
    push @foundtest, $v2; # I just want to test it.
    $vcs_found_aref = \@foundtest;
}
foreach my $e (@foundtest) {
    my $line = $e->chr . "\t" . $e->ref . "\t" . $e->pos . "\t" . $e->alt;
    print "FOUND VCS: $line\n";
}
# END TEST

# Step 3: Compute strand support
my $not_found_results_aref = &computeStrandSupport($vcs_not_found_aref, $pileup_fpath, "not found vcs");
my $found_results_aref = &computeStrandSupport($vcs_found_aref, $pileup_fpath, "found vcs");


# Step 4: Store results into file
&writeToFile($not_found_results_aref, $not_found_outfile);
&writeToFile($found_results_aref, $found_outfile);

print "Completed: $not_found_outfile\n";
print "Completed: $found_outfile\n";

### SUB-ROUTINES ##########################################
sub writeToFile() {
    my ($results_aref, $outfile) = @_;
    print "Writing out to $outfile\n";

    my @results = @{$results_aref};
    open(OUT, ">$outfile") || die "Unable to open $outfile\n";
    my $header = join("\t", "#CHROM", "POS", "REF", "ALT", 
		      "CONSENSUS_QUAL", "SNP_QUALITY", "MAX_MAPPING_QUALITY", "NUM_READS_COVERING_SITE", 
		      "READS_SUPPORT_REF_PLUS", "READS_SUPPORT_REF_MINUS", 
		      "READS_SUPPORT_ALT_PLUS", "READS_SUPPORT_ALT_MINUS");
    print OUT "$header\n";
    foreach my $v (@results) {
	my $line = join("\t", $v->chr, $v->pos, $v->ref, $v->alt, 
			$v->consensus_qual,
			$v->snp_quality, 
			$v->max_mapping_quality, 
			$v->num_reads_covering_site,
			$v->reads_support_ref_plus,
			$v->reads_support_ref_minus,
			$v->reads_support_alt_plus,
			$v->reads_support_alt_minus);
	print "$line\n";
	print OUT "$line\n";
    }
    close(OUT);

}#end writeToFile




sub splitVCSFoundAndNotFound() {
    my ($vcs_found_href, $disease, $database, $table) = @_;

    my %found_vcs = %{$vcs_found_href};

    my @not_found_results = ();
    my @found_results = ();

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$database", "", "", { RaiseError => 1, AutoCommit => 1 } );
    $dbh->do('PRAGMA synchronous=1');
    $dbh->do('PRAGMA cache_size=4000');
    my $query = "SELECT CHROMOSOME, REF_ALLELE, GENOMIC_COORD, ALT_ALLELE FROM $table WHERE DISEASE=?";
    my $stm = $dbh->prepare($query);
    $stm->execute($disease);
    while(my @row_array = $stm->fetchrow_array()) {
	my ($chr, $ref, $pos, $alt) = @row_array;
	my $key = "$chr:$pos:$ref:$alt";

	my $v = VARIANT_STRAND_SUPPORT->new();
	if ($chr !~ /^chr/) { $chr = "chr" . $chr; }
	$v->chr($chr);
	$v->ref($ref);
	$v->pos($pos);
	$v->alt($alt);

	if (!defined($found_vcs{$key})) {
	    # Then this VCS was not found in the sample, we store it.
	    push @not_found_results, $v;
	} else {
	    push @found_results, $v;
	}
    }#end while(my @row_array = $stm->fetchrow_array())


    return (\@found_results, \@not_found_results);
}#end getVCSNotFound

sub getVCSFound() {
    my ($foundVariantsFile) = @_;

    my %results = ();
    
    open(FOUND, "<$foundVariantsFile") || die "Unable to open $foundVariantsFile\n";
    while(my $row = <FOUND>) {
	chomp $row;
	if ($row !~ /\#/) {
	    my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $IS_VCS, $GENE_NAME) = split(/\t/, $row);
	    if ($IS_VCS eq $VCS) { # Then this is a VCS
		if ($CHROM !~ /^chr/) { $CHROM = "chr" . $CHROM; }
		my $key = "$CHROM:$POS:$REF:$ALT";
		$results{$key} = $row;
	    }
	}
    }#end while
    close(FOUND);

    return \%results;

}#end getVCSFound



__END__





VARIANT_STRAND_SUPPORT
    chr => '$',
    transcript_id => '$',
    gene_name => '$',    
    phenotype => '$',
    g_coord_1 => '$',
    g_coord_2 => '$',
    mut_type => '$',
    ref_amino => '$',
    alt_amino => '$',
    ref_allele => '$',
    alt_allele => '$',
    ref_pos => '$',
    alt_pos => '$',

    # This we get from computing 
    consensus_qual => '$',
    snp_quality => '$', 
    max_mapping_quality => '$', 
    num_reads_covering_site => '$',
    reads_support_ref_plus => '$',
    reads_support_ref_minus => '$',
    reads_support_alt_plus => '$',
    reads_support_alt_minus => '$',








sub computeStrandSupport() {
    my ($vcs_aref, $pileup_fpath, $type) = @_;
    print "Computing strand support for $type\n";

    my @results = ();
    my @variants = @{$vcs_aref};

    foreach my $v (@variants) {
	my $chr = $v->chr;
	my $g_coord = $v->g_coord_1;
	my $ref = $v->ref_allele;
	my $alt = $v->alt_allele;
	my $ref_plus_strand = "."; my $ref_minus_strand = ",";
	my $alt_plus_strand = uc($alt); my $alt_minus_strand = lc($alt);

	my @query_result = `grep -P "$chr\t$g_coord\t" $pileup_fpath`;
	if (scalar(@query_result) > 0) {
	    my $row_of_interest = $query_result[0];
	    my ($chrom, $coord, $base1, $base2, 
		$consensus_qual, $snp_quality, $max_mapping_quality, $num_reads_covering_site,
		$sequence_of_reads, @others) = split(/\t/, $row_of_interest);
	    my ($ref_plus_cnt, $ref_minus_cnt, $alt_plus_cnt, $alt_minus_cnt) = 
		&parseReads($sequence_of_reads, $alt_plus_strand, $alt_minus_strand);
	    
	    $v->consensus_qual($consensus_qual);
	    $v->snp_quality($snp_quality);
	    $v->max_mapping_quality($max_mapping_quality);
	    $v->num_reads_covering_site($num_reads_covering_site);

	    $v->reads_support_ref_plus($ref_plus_cnt);
	    $v->reads_support_ref_minus($ref_minus_cnt);
	    $v->reads_support_alt_plus($alt_plus_cnt);
	    $v->reads_support_alt_minus($alt_minus_cnt);
	    print "$chr, $g_coord, $ref => [$ref_plus_cnt, $ref_minus_cnt], $alt => [$alt_plus_cnt, $alt_minus_cnt]\n";

	} else {
	    $v->reads_support_ref_plus(0);
	    $v->reads_support_ref_minus(0);
	    $v->reads_support_alt_plus(0);
	    $v->reads_support_alt_minus(0);
	    $v->consensus_qual(0);
	    $v->snp_quality(0);
	    $v->max_mapping_quality(0); 
	    $v->num_reads_covering_site(0);

	    print "Missing: $chr, $g_coord, $ref, $alt\n";
	    
	}#end if (scalar(@query_result) > 0)

    }#end foreach my $v (@variants)

    
    return \@results;
}#end computeStrandSupport

sub parseReads() {


=pod 
# from http://en.wikipedia.org/wiki/Pileup_format
# since samtools documentation is so sparse I can't find the description.

Column 5: The bases string [edit]
. (dot) means a base that matched the reference on the forward strand
, (comma) means a base that matched the reference on the reverse strand
AGTCN denotes a base that did not match the reference on the forward strand
agtcn denotes a base that did not match the reference on the reverse strand
+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases
-[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases
^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
$ (dollar) marks the end of a read segment
=cut



    my ($sequence_of_reads, $alt_plus_string, $alt_minus_string) = @_;
    my $ref_plus_count = ($sequence_of_reads =~ tr/\./\./);
    my $ref_minus_count = ($sequence_of_reads =~ tr/\,/\,/);
    my $alt_plus_count =  ($sequence_of_reads =~ tr/$alt_plus_string/$alt_plus_string/);
    my $alt_minus_count =  ($sequence_of_reads =~ tr/$alt_minus_string/$alt_minus_string/);

    return ($ref_plus_count, $ref_minus_count, $alt_plus_count, $alt_minus_count);
}#end parseReads
