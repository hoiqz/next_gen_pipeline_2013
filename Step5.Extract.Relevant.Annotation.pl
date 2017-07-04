#!/usr/bin/perl -w

###########################################################################
#
# We now extract the relevant information from VEP
# 
# @author: Hoi Qiangze & Sim Ngak Leng
# @date: 2013.10.28
#
###########################################################################

use strict;
use Cwd 'abs_path';
use File::Basename;
use Class::Struct;
use DBI;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.objects.pl";
require "$pwd/polaris.annotation.selector.pl";
require "$pwd/polaris.normal.tumor.utility.pl";


if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}
print "Script: $0\n";
my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;

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
my $insubdir = $subdirs{MARK_VAR_CALLERS};
my $indir = "$subrootdir/$insubdir";
my $infile = "$indir/$rootname.marked.CAG.vcs"; # This is a VCF v4.0 file
my $isVCF = $userConfig{VEP_INFILE_FORMAT};

my $infile_not_CAG = "$indir/$rootname.marked.not.CAG.vcs"; # This is a VCF v4.0 file


##### THIS PART IS TO CATER TO IAIN TAN'S TUMOR-NORMAL 
# Idea: Add 2 new keys to User Config file -> NORMAL_SAMPLE & NORMAL_PATH
# THREE SCENARIOS:
# (1) If the sample is a Normal sample, the NORMAL_SAMPLE value will be empty
# (2) If the run is a Tumor sample run BUT there is no Normal Sample to compare with, the NORMAL_SAMPLE will ALSO be empty
# (3) Otherwise, NORMAL_PATH points to the root directory of the normal sample variant call
#     and NORMAL_SAMPLE points to the name: eg. SGH.H226_S12.NSCLC
#     eg. /mnt/PolarisPool/Lung_Panel_2nd_run_/SGH.H226_S12.NSCLC
# We append mark.var.callers/$NORMAL_SAMPLE.marked.CAG.vcs
# and mark.var.callers/$NORMAL_SAMPLE.marked.not.CAG.vcs
#
# We rename $infile & $infile_not_CAG 
# (CAG means Clinically Actionable Genes, determines if the variant goes to Clinical Report or Discovery Report)
# We take them and and normal sample's CAG.vcs, and run VCFTools to get variants in this but NOT in Normal
# and write into $infile and $infile_not_CAG
# 
my $normalSamplePath = $userConfig{NORMAL_PATH};
my $normalSample = $userConfig{NORMAL_SAMPLE};
if (defined $normalSample && $normalSample ne "") {
    my $isec = $softwareTools{ISEC};
    &runFilterWithNormal($infile, $normalSamplePath, $normalSample, $isec, "CAG");     # in clinically actionable genes
    &runFilterWithNormal($infile, $normalSamplePath, $normalSample, $isec, "not.CAG"); # not in clinically actionable genes
}
##### END IAIN TAN'S TUMOR-NORMAL



# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{SELECTED_VEP_ANNOTATION};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }

my $outfile = "$outdir/$rootname.CAG.selected.annotations";
my $discardedfile = "$outdir/$rootname.CAG.discarded"; # These are discarded because they are, e.g. synonymous or intron (even if they have COSMIC id)


my $outfile_not_CAG = "$outdir/$rootname.not.CAG.selected.annotations";
my $discardedfile_not_CAG = "$outdir/$rootname.not.CAG.discarded"; # These are discarded because they are, e.g. synonymous or intron (even if they have COSMIC id)



# GET REFLEX VCS DATABASE
my $database = $references{REFLEX_VCS_DATABASE};
my $db_table = $references{REFLEX_VCS_DATABASE_TABLE};
my $transcripts_href = &getTranscriptsForTestPanel($database, $db_table, $disease);


&extractForPOLARIS($discardedfile, $outfile, $infile, $transcripts_href);

&extractForPOLARIS($discardedfile_not_CAG, $outfile_not_CAG, $infile_not_CAG, $transcripts_href);

if ($DEBUG == $TRUE) {
    `cat /home/polaris/t1/VIR/injection/INDEL.TEST.CAG.selected.annotation >> $outfile`;
    `cat /home/polaris/t1/VIR/injection/INDEL.TEST.not.CAG.selected.annotation >> $outfile_not_CAG`;
}






sub getTranscriptsForTestPanel() {
    my ($database, $table, $disease) = @_;
    my %results = ();
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$database", "", "", { RaiseError => 1, AutoCommit => 1 } );
    $dbh->do('PRAGMA synchronous=1');
    $dbh->do('PRAGMA cache_size=4000');
    my $query = "SELECT DISTINCT TRANSCRIPT_ID FROM $table WHERE DISEASE=?";
    my $stm = $dbh->prepare($query);
    $stm->execute($disease);
    while(my @row_array = $stm->fetchrow_array()) {
	my $transcript = $row_array[0];
	$transcript =~ s/^\s+//;
	$transcript =~ s/\s+$//;
	$results{$transcript} = 0;
    }#end while(my @row_array = $stm->fetchrow_array())

    $dbh->disconnect;

    return \%results;
}# getTranscriptsForTestPanel



sub extractForPOLARIS() {

    my ($discardedfile, $outfile, $infile, $vcs_transcripts_href) = @_;

    open(DIS, ">$discardedfile") || die "Unable to open $discardedfile for writing.\n";
    open(OUT, ">$outfile") || die "Unable to open $outfile for writing.\n";
    open(IN, "<$infile") || die "Unable to open $infile for reading.\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row =~ /^\#/) {
	    print ">>> $row\n";
	    print OUT "$row\n";
	    print DIS "$row\n";
	}#end if ($row !~ /^\#/)
	else {
#	print "ROW: $row\n";
	    my @data = split(/\t/, $row);	
	    
	    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $value, $is_vcs, $gene_name) = split(/\t/, $row);
	    if ($isVCF !~ /vcf/i) {
		print "Is VarScan\n";
		($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $is_vcs, $gene_name) = split(/\t/, $row);
	    } else {
		print "Is MiSeq VCF\n";
	    }
	    
	    
	    
	    if ($info =~ /(.*)CSQ\=(.*)/) {
		my $others = $1;
		my $line = $2;
		# print "$info\n";
		# C|ENSG00000171094|ENST00000389048|Transcript|missense_variant|5288|4381|1461|I/V|Atc/Gtc|rs1670283|NM_004304.4|29/29||||T:0.0069|CCDS33172.1||0.03|0|0|0,
		# C|ENSG00000171094|ENST00000431873|Transcript|missense_variant|871 |871 |291 |I/V|Atc/Gtc|rs1670283|           |4/4  ||||T:0.0069|           ||0.03|0|0|0
		my @data = split(/,/, $line);
		
		my @candidates = ();
		foreach my $can (@data) { push @candidates, &storeAsObject($can); }
		my $selected = &selectCandidates(\@candidates, $vcs_transcripts_href);

		my $selected_annotation = "CSQ=" . &deflateIntoLine($selected);

		print "SELECTED ANNOTATION: $selected_annotation\n";

		if ($others && $others ne "") { $selected_annotation = "$others;$selected_annotation"; }
#	    print "IS_VCS: $is_vcs, GENE: $gene_name\n";
		my $new_row = join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $selected_annotation, $is_vcs, $gene_name);
		# Perform filtering
		# Filter out if 
		# not VCS && selected annotation has no COSMIC && selected annotation is not LoF or non-synonymous or not in exonic region
		my $filterOut = &filterVariants($selected->consequences, $is_vcs);
		if ($filterOut == $TRUE) {
		    print OUT "$new_row\n";
#		print "FILTERED: $new_row\n";
		} else {
		    print DIS "$row\n"; # I print out the original row just to be sure, easier to debug in the future.
		}
	    } else {
		# There's no VEP annotation at all, we throw it out.
		print DIS "$row\n";
	    }
	}#end else
    }#end while
    close(IN);
    close(OUT);
    close(DIS);
    
    print "Completed: $outfile\n";

} #end extractForPOLARIS()


### SUB-ROUTINES

sub storeAsObject() {
    my $data = $_[0];
    my @elements = split(/\|/, $data);
    #print "Number: " . scalar(@elements) . "\n";

    my $obj = VEP_ANNO->new();
    $obj->allele($elements[0]);
    $obj->gene_id($elements[1]);
    $obj->transcript_id($elements[2]);
    $obj->variant_type($elements[3]);
    $obj->consequences($elements[4]);
    $obj->cDNA_pos($elements[5]);
    $obj->cds_pos($elements[6]);
    $obj->protein_pos($elements[7]);
    $obj->amino_acids($elements[8]);
    $obj->codon_change($elements[9]);
    $obj->existing_variation($elements[10]);
    $obj->refseq($elements[11]);
    $obj->exon($elements[12]);
    $obj->intron($elements[13]);
    $obj->distance($elements[14]);
    $obj->clin_sig($elements[15]);
    $obj->gmaf($elements[16]);
    $obj->ccds($elements[17]);
    $obj->cell_type($elements[18]);
    $obj->afr_maf($elements[19]);
    $obj->amr_maf($elements[20]);
    $obj->asn_maf($elements[21]);
    $obj->eur_maf($elements[22]);
    $obj->priority_score(0); # Start off with no score
    return $obj;
}#end storeAsObject

sub deflateIntoLine() {
    # Takes a VEP_ANNO object and turn it into a VEP CSQ line
    my $obj = $_[0];


    if (!$obj->allele || $obj->allele eq "") { $obj->allele("."); }
    if (!$obj->gene_id || $obj->gene_id eq "") { $obj->gene_id("."); }
    if (!$obj->transcript_id || $obj->transcript_id eq "") { $obj->transcript_id("."); }
    if (!$obj->variant_type || $obj->variant_type eq "") { $obj->variant_type("."); }
    if (!$obj->consequences || $obj->consequences eq "") { $obj->consequences("."); }
    if (!$obj->cDNA_pos || $obj->cDNA_pos eq "") { $obj->cDNA_pos("."); }
    if (!$obj->cds_pos || $obj->cds_pos eq "") { $obj->cds_pos("."); }
    if (!$obj->protein_pos || $obj->protein_pos eq "") { $obj->protein_pos("."); }
    if (!$obj->amino_acids || $obj->amino_acids eq "") { $obj->amino_acids("."); }
    if (!$obj->codon_change || $obj->codon_change eq "") { $obj->codon_change("."); }
    if (!$obj->existing_variation || $obj->existing_variation eq "") { $obj->existing_variation("."); }
    if (!$obj->refseq || $obj->refseq eq "") { $obj->refseq("."); }
    if (!$obj->exon || $obj->exon eq "") { $obj->exon("."); }
    if (!$obj->intron || $obj->intron eq "") { $obj->intron("."); }
    if (!$obj->distance || $obj->distance eq "") { $obj->distance("."); }
    if (!$obj->clin_sig || $obj->clin_sig eq "") { $obj->clin_sig("."); }
    if (!$obj->gmaf || $obj->gmaf eq "") { $obj->gmaf("."); }
    if (!$obj->ccds || $obj->ccds eq "") { $obj->ccds("."); }
    if (!$obj->cell_type || $obj->cell_type eq "") { $obj->cell_type("."); }
    if (!$obj->afr_maf || $obj->afr_maf eq "") { $obj->afr_maf("."); }
    if (!$obj->amr_maf || $obj->amr_maf eq "") { $obj->amr_maf("."); }
    if (!$obj->asn_maf || $obj->asn_maf eq "") { $obj->asn_maf("."); }
    if (!$obj->eur_maf || $obj->eur_maf eq "") { $obj->eur_maf("."); }

    my $line = join("\|", 
		    $obj->allele,
		    $obj->gene_id,
		    $obj->transcript_id,
		    $obj->variant_type,
		    $obj->consequences,
		    $obj->cDNA_pos,
		    $obj->cds_pos,
		    $obj->protein_pos,
		    $obj->amino_acids,
		    $obj->codon_change,
		    $obj->existing_variation,
		    $obj->refseq,
		    $obj->exon,
		    $obj->intron,
		    $obj->distance,
		    $obj->clin_sig,
		    $obj->gmaf,
		    $obj->ccds,
		    $obj->cell_type,
		    $obj->afr_maf,
		    $obj->amr_maf,
		    $obj->asn_maf,
		    $obj->eur_maf);

    return $line;
    
}#end deflateIntoLine



__END__
