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
use Class::Struct;
use DBI;
#use XML::Writer;

my $absPath = abs_path($0);
my $pwd = dirname($absPath);
require "$pwd/readConfig.pl";
require "$pwd/polaris.varscan.support.pl";
require "$pwd/polaris.clinical.text.printer.pl";

if (scalar(@ARGV) != 4) {
    die "Usage: perl $0 <software config file> <subdir config file> <references file> <user config file>\n" .
	"Eg. perl $0 config/POLARIS.SOFTWARE.CONF config/polaris.subdirs.conf config/POLARIS.REFERENCE.CONF user.config\n";
}
print "Script: $0\n";
my ($swconfigfile, $subdirconfigfile, $referencesfile, $userconfigfile) = @ARGV;

my $TRUE = 0;
my $FALSE = 1;
my $NEED_JOURNAL_REFERENCES = $TRUE;

struct MUTATION => {
    # What is to be printed out
    gene_name => '$',
    alteration => '$', # e.g. G1244D, indel
    mutation_frequency => '$',
    annotation => '$', # This is dependent on whether the variant is a VCS or not.

    # These will be appended together with Annotation if they exist
    cosmic_id => '$',
    confirmed_somatic => '$', # I need to get this from COSMIC
    polyphen => '$',
    dbsnp => '$',
    esp => '$',
    g1000 => '$',

    # Internal structure
    is_validated => '$',
    is_vcs => '$',
    chr => '$',
    ref => '$',
    alt => '$',
    pos => '$',
    consequence => '$',

    # Other information from VarScan (if exists)
    # These information will not exist if the data comes from MiSeq's VCF
    Reads1 => '$', 
    Reads2 => '$', 
    VarFreq => '$', # This is the same as mutation_frequency, but I'll just add it in
    Strands1 => '$', 
    Strands2 => '$', 
    Qual1 => '$', 
    Qual2 => '$', 
    Pvalue => '$', 
    MapQual1 => '$', 
    MapQual2 => '$', 
    Reads1Plus => '$', 
    Reads1Minus => '$', 
    Reads2Plus => '$', 
    Reads2Minus => '$', 
    VarAllele => '$',
}; #end struct

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

my $is_vcf = $userConfig{VEP_INFILE_FORMAT};

# INPUT INFORMATION
my $insubdir = $subdirs{CHECK_INDEPENDENT_VALIDATION_DIR};
my $indir = "$subrootdir/$insubdir";
my $infile = "$indir/$rootname.indp.validate.not.CAG.annotated"; # This is the not Clinically Actionable Gene list
# This is a temporary workaround until Clinicians gives us the actual interpretation data
my $generic_interpretation_file = $references{GENERIC_GENE_DESCRIPTION_FILE};

my $diagnosis = $userConfig{DIAGNOSIS};
if (! $diagnosis) { $diagnosis = $userConfig{DISEASE}; }

# DATE
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
#my $date = strftime "%d %m %Y", localtime;
$year += 1900;
my $RUN_DATE = "$mday $months[$mon] $year";


# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{REFLEX_DISCOVERY_REPORTS};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $reflex_discovery_xml = "$outdir/$rootname.discovery.report.xml";
my $reflex_discovery_txt = "$outdir/$rootname.discovery.report.txt";
my $reflex_pdf = "$outdir/$rootname.discovery.pdf";




# Step 1: Collate all the information required (gene name, gene alteration, mutation frequency, annotations, literature references)
my $var_for_printing_aref;
if ($is_vcf =~ /vcf/i) {
    my $limited_dir = $subdirs{LIMITED_DIR};
    my $original_vcf = "$subrootdir/$limited_dir/$rootname.not.in.clinically.actionable.genes";
    $var_for_printing_aref = &generateFromVCF($infile, $original_vcf); # This is used when the variants come from MiSeq VCF
} else { # The variants were called using samtools pileup
    # I need the original Varscan
    my $varscansubdir = $subdirs{VARSCAN_DIR};
    my $varscandir = "$subrootdir/$varscansubdir";
    my $varscanfile = "$varscandir/$rootname.not.CAG.variants";
    $var_for_printing_aref = &generateFromVarScan($infile, $varscanfile);    
}#end handling VarScan-based variant calling


my $var_to_print_out_aref = &generateMoreInfoForDiscovery($var_for_printing_aref);


my %for_printing = ();

my $generic_database = $references{POLARIS_GENERIC_INTERPRETATION_DB};
my $generic_table = $references{POLARIS_GENERIC_INTERPRETATION_TABLE};

my ($num_variants, $VARIANTS_AREF) = &generateTableList($var_to_print_out_aref);

$for_printing{SAMPLE_ID_1} = $hospital;
$for_printing{SAMPLE_ID_2} = $patientid;
$for_printing{VARIANTS} = $VARIANTS_AREF;
$for_printing{NUM_VARIANTS} = $num_variants;


&printText(\%for_printing, $reflex_discovery_txt);

#&printXML(\%for_printing, $reflex_discovery_xml, $num_variants);
&printDiscoveryXML(\%for_printing, $reflex_discovery_xml, $num_variants);
if (-e $reflex_discovery_xml) {
    print "Completed: $reflex_discovery_xml\n";
} else {
    print "ERROR: XML file not generated: $reflex_discovery_xml\n";
}

print "Completed: $reflex_discovery_txt\n";



### THIS IS THE FINAL STEP, WHERE THE PDF FILE IS GENERATED
my $java_exec = $softwareTools{JAVA_EXEC};
my $jarfile = $references{REFLEX_PDF_JAR};
my $discovery_root_dir = $references{REFLEX_DISCOVERY_DIR};
my $jasper_file = $discovery_root_dir . "/" .  $references{REFLEX_DISCOVERY_JASPER};
my $jrxml_file = $discovery_root_dir . "/" .  $references{REFLEX_DISCOVERY_JRXML};
my $jrprint_file = $discovery_root_dir . "/" .  $references{REFLEX_DISCOVERY_JRPRINT};
my $jasper_subreportdir = $discovery_root_dir . "/subreport/";
my $pdf_output_file = "$outdir/$rootname.discovery.report.pdf";
my $pdf_cmd = "$java_exec -jar $jarfile -x $reflex_discovery_xml -j $jasper_file -p $jrprint_file -r $jrxml_file -s $jasper_subreportdir -o $pdf_output_file";
print "$pdf_cmd\n";
my $result = system($pdf_cmd);
if ($result == 0) {
    print "Succeeded in pdf creation: $pdf_output_file\n";
} else {
    print "Failed creation of pdf file\n";
    die "$pdf_cmd\n";
}


























### THIS IS FOR HUMAN READABILITY





### SUB-ROUTINES 
sub printText() {
    print "printText\n";
    my ($data_href, $outfile) = @_;
    my %data = %{$data_href};
    
    my $sample_id_1 = $data{SAMPLE_ID_1};
    my $sample_id_2 = $data{SAMPLE_ID_2};

    my $num_variants = $data{NUM_VARIANTS};

    my $header = join("\t", "SAMPLE_ID_1", "SAMPLE_ID_2", "IS_INDEPENDENTLY_VALIDATED",
		      "GENE", "MUTATION", "MUTATION_FREQUENCY", "COSMIC", "DBSNP", "ESP", "POLYPHEN", "1000 GENOME", "CONFIRMED_SOMATIC");

    my $variants_aref = $data{VARIANTS};
    my @variants = @{$variants_aref};
    
    open(OUT, ">$outfile") || die "Unable to open $outfile\n";
    print OUT "$header\n";
    if ($num_variants > 0) {
	foreach my $v (@variants) {
	    my @row = split(/\t/, $v);
	    my $info_for_txt_file = $row[-1];
	    my ($gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen) = split(/:/, $info_for_txt_file);
	    my $line = join("\t", $sample_id_1, $sample_id_2, "Yes", $gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen);
	    print OUT "$line\n";
	}#end foreach
    }

    close(OUT);
}#end printText


# To be superceded by printDiscoveryXML
sub printXML() {
    print "printXML - to change to using XML::Writer\n";
    my ($data_href, $outfile, $num_variants) = @_;
    my %data = %{$data_href};
    
    my $version = &getPipelineVersion();

    open(OUT, ">$outfile") || die "Unable to open $outfile for writing.\n";
    print OUT "<?xml version=\"1.0\" ?>\n";
    print OUT "<POLARIS>\n";
    print OUT "<BIOINFO_VERSION>$version</BIOINFO_VERSION>\n";
    print OUT "<RUN_DATE>$RUN_DATE</RUN_DATE>\n";
    print OUT "<DIAGNOSIS>$diagnosis</DIAGNOSIS>\n";
    print OUT "<SAMPLE_ID_1>$hospital</SAMPLE_ID_1>\n";
    print OUT "<SAMPLE_ID_2>$patientid</SAMPLE_ID_2>\n";
    
    # Summary line
#    my $summary_line = $data{SUMMARY};
#    print OUT "<SUMMARY>$summary_line</SUMMARY>\n";

    # Table 1
    print OUT "<GENOMIC_ALTERATIONS_FOUND>\n";

    if ($num_variants > 0) {

	my $table1_aref = $data{VARIANTS};
	my @table1 = @{$table1_aref};
	
	foreach my $mutation (@table1) {
	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, 
		$hasPolyPhen, $hasCOSMIC, $hasDBSNP, $hasESP, $has1000Genome, $confirmed_somatic, @others) = split(/\t/, $mutation);
	    print OUT "<GENE_ALTERATION>$gene_alteration</GENE_ALTERATION>\n";
	    
	    print OUT "<POLYPHEN>$hasPolyPhen</POLYPHEN>\n";
	    print OUT "<COSMIC>$hasCOSMIC</COSMIC>\n";
	    print OUT "<CONFIRMED_SOMATIC>$confirmed_somatic</CONFIRMED_SOMATIC>\n";
	    print OUT "<DBSNP>$hasDBSNP</DBSNP>\n";
	    print OUT "<ESP>$hasESP</ESP>\n";
	    print OUT "<GENOME1000>$has1000Genome</GENOME1000>\n";
	    print OUT "<MUTATION_FREQUENCY>$mut_freq</MUTATION_FREQUENCY>\n";
	    print OUT "<ANNOTATION>$annotation\n";
	    if ($gene_annotations ne "NONE") {
		print OUT "$gene_annotations\n";
	    }
	    print OUT "</ANNOTATION>\n";
	}

    } else {
	# We do not have anything to print out
	print OUT "<GENE_ALTERATION> </GENE_ALTERATION>\n";
	print OUT "<POLYPHEN> </POLYPHEN>\n";
	print OUT "<COSMIC> </COSMIC>\n";
	print OUT "<CONFIRMED_SOMATIC> </CONFIRMED_SOMATIC>\n";
	print OUT "<DBSNP> </DBSNP>\n";
	print OUT "<ESP> </ESP>\n";
	print OUT "<GENOME1000> </GENOME1000>\n";
	print OUT "<MUTATION_FREQUENCY> </MUTATION_FREQUENCY>\n";
	print OUT "<ANNOTATION>NA</ANNOTATION>\n";


    }
    print OUT "</GENOMIC_ALTERATIONS_FOUND>\n";
    
    
    print OUT "</POLARIS>\n";
    close(OUT);
}#end printXML


sub printDiscoveryXML() {
    print "printDiscoveryXML\n";
    my ($data_href, $outfile, $num_variants) = @_;
    my %data = %{$data_href};

    my $version = &getPipelineVersion();

    open(OUT, ">$outfile") || die "Unable to open $outfile for writing.\n";
    print OUT "<?xml version=\"1.0\" ?>\n";
    print OUT "<POLARIS>\n";
    print OUT "<BIOINFO_VERSION>$version</BIOINFO_VERSION>\n";
    print OUT "<RUN_DATE>$RUN_DATE</RUN_DATE>\n";
    print OUT "<DIAGNOSIS>$diagnosis</DIAGNOSIS>\n";
    print OUT "<SAMPLE_ID_1>$hospital</SAMPLE_ID_1>\n";
    print OUT "<SAMPLE_ID_2>$patientid</SAMPLE_ID_2>\n";
    
    # Summary line
#    my $summary_line = $data{SUMMARY};
#    print OUT "<SUMMARY>$summary_line</SUMMARY>\n";
    print OUT "<SUMMARY> </SUMMARY>\n";

    # Table 1
    print OUT "<VALIDATED_GENOMIC_ALTERATIONS_FOUND>\n";
#    my $table1_aref = $data{IS_INDP_VALIDATED};
    my $table1_aref = $data{VARIANTS};

#    my $num_indp_validated = $data{NUM_INDP_VALIDATED};
#    if ($num_indp_validated == 0) {
    if ($num_variants == 0) {
#	print OUT "<GENE_ALTERATION>NONE</GENE_ALTERATION>\n";
	print OUT "<GENE_ALTERATION>NONE\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	print OUT "<MUTATION_FREQUENCY>NA</MUTATION_FREQUENCY>\n";
	print OUT "<ANNOTATION>NA</ANNOTATION>\n";	    
	print OUT "</GENE_ALTERATION>\n";
    } else { 
	my @table1 = @{$table1_aref};
	foreach my $mutation (@table1) {
###	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, @others) = split(/\t/, $mutation);
	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, 
		$hasPolyPhen, $hasCOSMIC, $hasDBSNP, $hasESP, $has1000Genome, $confirmed_somatic, @others) = split(/\t/, $mutation);

	    print "$mutation\n";
	    print "$has1000Genome\n";

	    print OUT "<GENE_ALTERATION>$gene_alteration\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	    print OUT "<MUTATION_FREQUENCY>$mut_freq</MUTATION_FREQUENCY>\n";
	    print OUT "<ANNOTATION>$annotation\n";
	    if ($gene_annotations ne "NONE") {
		if ($has1000Genome) {
		    print OUT "$gene_annotations,1000Genomes\n";
		} else {
		    print OUT "$gene_annotations\n";
		}
	    }
	    print OUT "</ANNOTATION>\n";
	    print OUT "</GENE_ALTERATION>\n";
	}
    }
    print OUT "</VALIDATED_GENOMIC_ALTERATIONS_FOUND>\n";


    # Table 2
    print OUT "<OTHER_GENOMIC_ALTERATIONS_FOUND>\n";
    my $table2_aref = $data{NOT_INDP_VALIDATED};
    my $num_not_indp_validated = $data{NUM_NOT_INDP_VALIDATED};
#    if ($num_not_indp_validated == 0) {
#	print OUT "<GENE_ALTERATION>NONE</GENE_ALTERATION>\n";
	print OUT "<GENE_ALTERATION>NONE\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	print OUT "<MUTATION_FREQUENCY> NA </MUTATION_FREQUENCY>\n";
	print OUT "<ANNOTATION>NA</ANNOTATION>\n";
	print OUT "</GENE_ALTERATION>\n";
#    } else {
#	my @table2 = @{$table2_aref};
#	foreach my $mutation (@table2) {
#	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, @others) = split(/\t/, $mutation);
#	    print OUT "<GENE_ALTERATION>$gene_alteration\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
#	    print OUT "<MUTATION_FREQUENCY>$mut_freq</MUTATION_FREQUENCY>\n";
#	    print OUT "<ANNOTATION>$annotation\n";
#	    if ($gene_annotations ne "NONE") {
#		print OUT "$gene_annotations\n";
#	    }
#	    print OUT "</ANNOTATION>\n";
#	    print OUT "</GENE_ALTERATION>\n";
#	}
#    }
    print OUT "</OTHER_GENOMIC_ALTERATIONS_FOUND>\n";

    # Now we print out methods
    my $methods = &getClinicalReportMethodsText($disease, $references_href);
    print OUT "<METHODS>\n";
    print OUT "$methods"; # DO NOT ADD NEWLINE, the variable represents more than 1 paragraph of text and already has newlines
    print OUT "</METHODS>\n";

    my $limitations = &getLimitationsText($disease, $references_href);
    print OUT "<LIMITATIONS>\n";
    print OUT "$limitations"; # DO NOT ADD NEWLINE, the variable represents more than 1 paragraph of text and already has newlines
    print OUT "</LIMITATIONS>\n";

    my $about = &getAboutText($disease, $references_href);
    print OUT "<ABOUT>\n";
    print OUT $about;
    print OUT "</ABOUT>\n";


    print OUT "</POLARIS>\n";
    close(OUT);

}#end printXML



































sub generateTableList() {
    my ($aref, $is_indp_validated) = @_;
    my @data = @{$aref};
    my @results = ();
    if (scalar(@data) == 0) {
	my $datum = join("\t", "NONE FOUND", "NA", "NA");	
	push @results, $datum;
	return (0, \@results);
    }

    foreach my $mut (@data) {

	my $gene = $mut->gene_name;
	my $alteration = $mut->alteration;

	if ($alteration eq "") {
	    # This is because Amino Acid change does not exist for splice, frameshifts, etc.
	    # So we just print out the consequence (frameshift, etc.)
	    my $conseq = $mut->consequence;
	    my @conseqs = split(/\&/, $conseq);
	    my $c = $conseqs[0];
	    my @variants = split(/\_/, $c);
	    my $joined = join(" ", @variants);
	    $alteration = $joined;
	}

	print "$gene ----> $alteration\n";
=pod
	if ($alteration eq "") {
	    # Check if it is splice
	    if ($mut->consequence =~ /splice/) {
		my $chr = $mut->chr; my $coord = $mut->pos; my $ref = $mut->ref; my $alt = $mut->alt;
		$alteration = "splice [$chr,$coord,$ref\/$alt]";
	    }
	}
=cut
	my $gene_alteration = "$gene $alteration";
	my $mut_freq = $mut->mutation_frequency;
	my $annotation = &getGenericInterpretation($gene, $generic_database, $generic_table);
	
	my $for_txt_file = join(":", $gene, $alteration, $mut_freq);
	
	my $confirmed_somatic = $mut->confirmed_somatic;
	if (! $confirmed_somatic) { $confirmed_somatic = "N"; }
	my $hasPolyPhen = "N";
	my $hasCOSMIC = "N";
	my $hasDBSNP = "N";
	my $hasESP = "N";
	my $has1000G = "N";
	
	
	my @gene_annotations = ();
	if ($mut->cosmic_id && $mut->cosmic_id =~ /COSM/) {
	    push @gene_annotations, "COSMIC: " . $mut->cosmic_id;
	    $for_txt_file .= ":" . $mut->cosmic_id;
	    $hasCOSMIC = "Y";
	} else { $for_txt_file .= ":NA"; }
	
	if ($mut->dbsnp && $mut->dbsnp =~ /rs/) {
	    push @gene_annotations, "dbSNP: " . $mut->dbsnp;
	    $for_txt_file .= ":" . $mut->dbsnp;
	    $hasDBSNP = "Y";
	} else { $for_txt_file .= ":NA"; }
	if ($mut->esp && $mut->esp =~ /ESP/) {
	    push @gene_annotations, "ESP: " . $mut->esp;
	    $for_txt_file .= ":" . $mut->esp;
	    $hasESP = "Y";
	} else { $for_txt_file .= ":NA"; }
	if ($mut->polyphen && $mut->polyphen ne "." && $mut->polyphen ne "") {
	    push @gene_annotations, "PolyPhen: " . $mut->polyphen;
	    $for_txt_file .= ":" . $mut->polyphen;
	    $hasPolyPhen = "Y";
	} else { $for_txt_file .= ":NA"; }
	
	my $has1000Genome = $mut->g1000;
	
	
	my $gene_anno_line = "";
	if (scalar(@gene_annotations) == 0) {
	    $gene_anno_line = "NONE";
	} else {
	    $gene_anno_line = join(", ", @gene_annotations);
	}
	
	my $datum = join("\t", $gene_alteration, $mut_freq, $annotation, $gene_anno_line, $hasPolyPhen, $hasCOSMIC, $hasDBSNP, $hasESP, $has1000Genome, $confirmed_somatic, $for_txt_file);
	
	print "DATUM: $datum\n";
	push @results, $datum;
    }#end foreach

    my $num = scalar(@results);

    return ($num, \@results);
}#end

sub getGenericInterpretation() {
    print "getGenericInterpretation: This needs to be modified once disease champions provide\n" .
	"actual interpretations for non-validated variants in clinically actionable genes\n" .
	"Currently, we are using a tab-delimited file generated from RefSeq (ask Hoi how it was done).\n";
    my ($gene, $generic_database, $generic_table) = @_;

    my $gene_col = '$1'; my $interpret_col = '$13';
    my @data = `awk -F"\t" '{ if ($gene_col ~ /$gene/) { print $interpret_col }}' $generic_database`;
    my $interpretation = "";
    if (scalar(@data) > 0) {	
	$interpretation = $data[0];
    }
    if (!$interpretation || $interpretation eq "") { $interpretation = "No annotations available."; }
    chomp $interpretation;
    return $interpretation;
}#end getGenericInterpretation


sub generateSummaryLine() {
    my ($aref) = @_;
    my @data = @{$aref};

    # Handle case where there are no variants found in clinically actionable genes
    if (scalar(@data) == 0) {
	return "NO MUTATIONS FOUND IN CLINICALLY ACTIONABLE GENES";
    }
    
    # Otherwise, there is information to be printed out
    my @result = ();
    foreach my $mutation (@data) {
	my $gene = $mutation->gene_name;
	my $alteration = $mutation->alteration;
	push @result, "$gene $alteration";
    }

    my $line = join(",", @result) . " MUTATIONS";
    return $line;
}#end generateSummaryLine







### SUB-ROUTINES HANDLING MISEQ VCF ######################################################################
sub generateFromVCF() {
    print "generateFromVCF --> to be implemented 2013.10.30\n";
    my ($infile, $orig_vcf_file) = @_;
    my @results = ();


    # What I have here is a VCF file
    # First, get variant frequencies from original VCF file
    my $variant_frequencies_href = &parseVCF_to_get_VarFreq($orig_vcf_file);
    my %variant_frequencies = %{$variant_frequencies_href};

    open(IN, "<$infile") || die "Unable to open $infile for reading.\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row !~ /\#/) {

	    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $is_vcs, $gene_name, $indp_validated) = split(/\t/, $row);
	    my $key = "$chr:$pos:$ref:$alt";
	    my $vf = $variant_frequencies{$key};
	    if (!defined($vf)) { $vf = "0.00%"; }
	    my $mut = &createVCFMutation($row, $vf); # This is the MiSeq equivalent of createVarScanMutation
	    push @results, $mut;
	    print "MUTATION FREQ: $key - " . $mut->mutation_frequency . "\n";
	    print "GENE NAME: $key - " . $mut->gene_name . "\n";

	}#end if ($row !~ /\#/)
    }#end while
    close(IN);
    
    return \@results;
}#end generateFromVCF


sub createVCFMutation() {
    my ($row, $vf) = @_;

    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $is_vcs, $gene_name, $indp_validated) = split(/\t/, $row);

    # At this point, we have what we need to populate the mutation object
    my $mut = MUTATION->new();
    if ($chr !~ /^chr/) { $chr = "chr" . $chr; }
    $mut->chr($chr);
    $mut->gene_name($gene_name);
    print "GeneName: $gene_name\n";
    $mut->ref($ref);
    $mut->alt($alt);
    $mut->pos($pos);

    $info =~ /.*\;CSQ\=(.*)/;
    my $csq = $1;
    my ($Allele, $Gene, $Feature, $Feature_type, $Consequence, $cDNA_position, $CDS_position, 
	$Protein_position, $Amino_acids, $Codons, $Existing_variation, $RefSeq, 
	$EXON, $INTRON, $DISTANCE, $CLIN_SIG, $PolyPhen, 
	$GMAF, $CCDS, $CELL_TYPE, $AFR_MAF, $AMR_MAF, $ASN_MAF, $EUR_MAF) = split(/\|/, $csq);

    print "from Step8.Generate.Discovery.Report.pl\n";
    my $alteration = &generateAlteration($Protein_position, $Amino_acids, $ref, $alt, $chr, $csq);

    $mut->alteration($alteration); # e.g. G1244D, indel
    $mut->mutation_frequency($vf);

    $mut->annotation("."); # This is dependent on whether the variant is a VCS or not.

    my $cosmic_id = &getAnnotations($Existing_variation, "COSM");
    my $rsid = &getAnnotations($Existing_variation, "rs");
    my $esp = &getAnnotations($Existing_variation, "ESP");
    $mut->cosmic_id($cosmic_id);
    $mut->polyphen($PolyPhen);
    $mut->dbsnp($rsid);
    $mut->esp($esp);

    $mut->is_validated($indp_validated);
    $mut->is_vcs($is_vcs);

=pod MiSeq VCF does not provide us with these information
    $mut->Reads1($Reads1); 
    $mut->Reads2($Reads2); 
    $mut->VarFreq($VarFreq); # This is the same as mutation_frequency, but I'll just add it in
    $mut->Strands1($Strands1); 
    $mut->Strands2($Strands2); 
    $mut->Qual1($Qual1); 
    $mut->Qual2($Qual2); 
    $mut->Pvalue($Pvalue); 
    $mut->MapQual1($MapQual1); 
    $mut->MapQual2($MapQual2); 
    $mut->Reads1Plus($Reads1Plus); 
    $mut->Reads1Minus($Reads1Minus); 
    $mut->Reads2Plus($Reads2Plus); 
    $mut->Reads2Minus($Reads2Minus); 
    $mut->VarAllele($VarAllele);
=cut

    return $mut;

}#end createVCFMutation


sub parseVCF_to_get_VarFreq() {
    my ($infile) = @_;
    my %results = ();
    open(IN, "<$infile") || die "Unable to open $infile for reading.\n";
    while(my $row = <IN>) {
	chomp $row;
	if ($row !~ /\#/) {
	    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $values) = split(/\t/, $row);
	    my $key = "$chr:$pos:$ref:$alt";
	    my $variant_frequency = 0.00; # this should not be possible
	    my @formats = split(/:/, $format);
	    my $index = 0; # I need to search for the position of VF (cannot assume it's always the last
	    for(my $i = 0; $i < scalar(@formats); $i++) {
		my $f = $formats[$i];
		if ($f eq "VF") {
		    $index = $i;
		    last;
		}
	    }
	    my @vals = split(/:/, $values);
	    if ($index > 0) { $variant_frequency = $vals[$index]; }
	    my $var_freq = $variant_frequency * 100.0;
	    $results{$key} = "$var_freq" . "%";
	}#end if ($row !~ /\#/)
    }#end while
    close(IN);
    return \%results;
}#end parseVCF_to_get_VarFreq



__END__





