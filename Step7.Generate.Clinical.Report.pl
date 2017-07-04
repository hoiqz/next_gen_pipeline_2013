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
use XML::Writer;

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
my $infile = "$indir/$rootname.indp.validate.CAG.annotated";
# This is a temporary workaround until Clinicians gives us the actual interpretation data
my $generic_interpretation_file = $references{GENERIC_GENE_DESCRIPTION_FILE};


my $diagnosis = $userConfig{DIAGNOSIS};
if (! $diagnosis || $diagnosis eq "") { $diagnosis = $userConfig{DISEASE}; }

# DATE
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
#my $date = strftime "%d %m %Y", localtime;
$year += 1900;
my $RUN_DATE = "$mday $months[$mon] $year";




# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{REFLEX_REPORTS};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) { system("mkdir -p $outdir"); }
my $reflex_validated_xml = "$outdir/$rootname.clinical.report.xml";
my $reflex_pdf = "$outdir/$rootname.pdf";

# Step 1: Collate all the information required (gene name, gene alteration, mutation frequency, annotations, literature references)
my $var_to_print_out_aref;
if ($is_vcf =~ /vcf/i) {
    my $limited_dir = $subdirs{LIMITED_DIR};
    my $original_vcf = "$subrootdir/$limited_dir/$rootname.limited.to.clinically.actionable.genes";
    $var_to_print_out_aref = &generateFromVCF($infile, $original_vcf); # This is used when the variants come from MiSeq VCF
} else { # The variants were called using samtools pileup
    # I need the original Varscan
    my $varscansubdir = $subdirs{VARSCAN_DIR};
    my $varscandir = "$subrootdir/$varscansubdir";
    my $varscanfile = "$varscandir/$rootname.CAG.variants";
    $var_to_print_out_aref = &generateFromVarScan($infile, $varscanfile);
}#end handling VarScan-based variant calling

# Now generate lines to be used
my $found_mutations = &generateSummaryLine($var_to_print_out_aref);
print "FOUND MUTATIONS $found_mutations\n";


my %for_printing = ();
$for_printing{SUMMARY} = $found_mutations;

my $vcs_database = $references{POLARIS_VALIDATED_DB};
my $vcs_table = $references{POLARIS_VALIDATED_TABLE};

my $generic_database = $references{POLARIS_GENERIC_INTERPRETATION_DB};
my $generic_table = $references{POLARIS_GENERIC_INTERPRETATION_TABLE};


my $dbh = DBI->connect( "dbi:SQLite:dbname=$vcs_database", "", "", { RaiseError => 1, AutoCommit => 1 } );
$dbh->do('PRAGMA synchronous=1');
$dbh->do('PRAGMA cache_size=4000');
my $query = "SELECT INTERPRETATION FROM $vcs_table WHERE DISEASE=? AND CHROMOSOME=? AND REF_ALLELE=? AND GENOMIC_COORD=? AND ALT_ALLELE=?";
my $stm = $dbh->prepare($query);


my ($num_indp_validated, $IS_INDP_VALIDATED_TABLE_AREF) = &generateTableList($var_to_print_out_aref, "IS_INDP_VALIDATED");
my ($num_not_indp_validated, $NOT_INDP_VALIDATED_TABLE_AREF) = &generateTableList($var_to_print_out_aref, "NOT_INDP_VALIDATED");


$dbh->disconnect();

$for_printing{SAMPLE_ID_1} = $hospital;
$for_printing{SAMPLE_ID_2} = $patientid;
$for_printing{IS_INDP_VALIDATED} = $IS_INDP_VALIDATED_TABLE_AREF;
$for_printing{NOT_INDP_VALIDATED} = $NOT_INDP_VALIDATED_TABLE_AREF;
$for_printing{NUM_INDP_VALIDATED} = $num_indp_validated;
$for_printing{NUM_NOT_INDP_VALIDATED} = $num_not_indp_validated;


my $journal_dir=$subdirs{MARK_VAR_CALLERS};

my $journal_references = "$subrootdir/$journal_dir/$rootname.journal.references.CAG.txt";


&printXML(\%for_printing, $reflex_validated_xml, $journal_references);
if (-e $reflex_validated_xml) {
    print "Completed: $reflex_validated_xml\n";
} else {
    print "ERROR: XML file not generated: $reflex_validated_xml\n";
}


=pod
# TESTING
my @test1 = @{$IS_INDP_VALIDATED_TABLE_AREF};
foreach my $line (@test1) {
    print "IV: $line\n";
}
my @test2 = @{$NOT_INDP_VALIDATED_TABLE_AREF};
foreach my $line (@test2) {
    print "NOT: $line\n";
}#end foreach
# END TESTING
=cut



### THIS IS FOR HUMAN READABILITY
my $reflex_validated_txt = "$outdir/$rootname.clinical.report.txt";
&printText(\%for_printing, $reflex_validated_txt);
print "Completed: $reflex_validated_txt\n";

### THIS IS THE FINAL STEP, WHERE THE PDF FILE IS GENERATED
my $java_exec = $softwareTools{JAVA_EXEC};
my $jarfile = $references{REFLEX_PDF_JAR};
my $clinical_root_dir = $references{REFLEX_CLINICAL_DIR};
my $jasper_file = $clinical_root_dir . "/" .  $references{REFLEX_CLINICAL_JASPER};
my $jrxml_file = $clinical_root_dir . "/" .  $references{REFLEX_CLINICAL_JRXML};
my $jrprint_file = $clinical_root_dir . "/" .  $references{REFLEX_CLINICAL_JRPRINT};
my $jasper_subreportdir = $clinical_root_dir . "/subreport/";
my $pdf_output_file = "$outdir/$rootname.clinical.report.pdf";
my $pdf_cmd = "$java_exec -jar $jarfile -x $reflex_validated_xml -j $jasper_file -p $jrprint_file -r $jrxml_file -s $jasper_subreportdir -o $pdf_output_file";
print "$pdf_cmd\n";
my $result = system($pdf_cmd);
if ($result == 0) {
    print "Succeeded in pdf creation: $pdf_output_file\n";
} else {
    print "Failed creation of pdf file\n";
    die "$pdf_cmd\n";
}


### SUB-ROUTINES 
sub printText() {
    print "printText\n";
    my ($data_href, $outfile) = @_;
    my %data = %{$data_href};
    
    my $sample_id_1 = $data{SAMPLE_ID_1};
    my $sample_id_2 = $data{SAMPLE_ID_2};

    $num_indp_validated = $data{NUM_INDP_VALIDATED};
    $num_not_indp_validated = $data{NUM_NOT_INDP_VALIDATED};



    my $header = join("\t", "RUN_DATE", "SAMPLE_ID_1", "SAMPLE_ID_2", "IS_INDEPENDENTLY_VALIDATED",
		      "GENE", "MUTATION", "MUTATION_FREQUENCY", "COSMIC", "DBSNP", "ESP", "POLYPHEN");

    my $validated_aref = $data{IS_INDP_VALIDATED};
    my @validated = @{$validated_aref};
    my $unvalidated_aref = $data{NOT_INDP_VALIDATED};
    my @unvalidated = @{$unvalidated_aref};
    
    open(OUT, ">$outfile") || die "Unable to open $outfile\n";
    print OUT "$header\n";

    if ($num_indp_validated > 0) {
	foreach my $v (@validated) {
	    my @row = split(/\t/, $v);
	    my $info_for_txt_file = $row[-1];
	    print "$info_for_txt_file\n";
	    my ($gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen) = split(/:/, $info_for_txt_file);
	    my $line = join("\t", $RUN_DATE, $sample_id_1, $sample_id_2, "Yes", $gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen);
	    print OUT "$line\n";
	}#end foreach
    }
    if ($num_not_indp_validated > 0) {
	foreach my $v (@unvalidated) {
	    print "UNVALIDATED ";
	    my @row = split(/\t/, $v);
	    my $info_for_txt_file = $row[-1];
	    print "$info_for_txt_file\n";
	    my ($gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen) = split(/:/, $info_for_txt_file);
	    my $line = join("\t", $RUN_DATE, $sample_id_1, $sample_id_2, "No", $gene, $alteration, $mut_freq, $cosmic, $dbsnp, $esp, $polyphen);
	    print OUT "$line\n";
	}#end foreach
    }
    close(OUT);
}#end printText




sub printXML() {
    print "printXML - to change to using XML::Writer\n";
    my ($data_href, $outfile, $journal_references) = @_;
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
    my $summary_line = $data{SUMMARY};
    print OUT "<SUMMARY>$summary_line</SUMMARY>\n";

    # Table 1
    print OUT "<VALIDATED_GENOMIC_ALTERATIONS_FOUND>\n";
    my $table1_aref = $data{IS_INDP_VALIDATED};

    my $num_indp_validated = $data{NUM_INDP_VALIDATED};
    if ($num_indp_validated == 0) {
#	print OUT "<GENE_ALTERATION>NONE</GENE_ALTERATION>\n";
	print OUT "<GENE_ALTERATION>NONE\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	print OUT "<MUTATION_FREQUENCY>NA</MUTATION_FREQUENCY>\n";
	print OUT "<ANNOTATION>NA</ANNOTATION>\n";	    
	print OUT "</GENE_ALTERATION>\n";
    } else { 
	my @table1 = @{$table1_aref};
	foreach my $mutation (@table1) {
	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, @others) = split(/\t/, $mutation);
#	    print OUT "<GENE_ALTERATION>$gene_alteration</GENE_ALTERATION>\n";
	    print OUT "<GENE_ALTERATION>$gene_alteration\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	    print OUT "<MUTATION_FREQUENCY>$mut_freq</MUTATION_FREQUENCY>\n";
	    print OUT "<ANNOTATION>$annotation\n";
	    if ($gene_annotations ne "NONE") {
		print OUT "$gene_annotations\n";
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
    if ($num_not_indp_validated == 0) {
#	print OUT "<GENE_ALTERATION>NONE</GENE_ALTERATION>\n";
	print OUT "<GENE_ALTERATION>NONE\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	print OUT "<MUTATION_FREQUENCY> NA </MUTATION_FREQUENCY>\n";
	print OUT "<ANNOTATION>NA</ANNOTATION>\n";
	print OUT "</GENE_ALTERATION>\n";
    } else {
	my @table2 = @{$table2_aref};
	foreach my $mutation (@table2) {
	    my ($gene_alteration, $mut_freq, $annotation, $gene_annotations, @others) = split(/\t/, $mutation);
#	    print OUT "<GENE_ALTERATION>$gene_alteration</GENE_ALTERATION>\n"; 
	    print OUT "<GENE_ALTERATION>$gene_alteration\n"; # 2013.12.02 To accommodate Seow Eng's bug fix of only 1 row showing up.
	    print OUT "<MUTATION_FREQUENCY>$mut_freq</MUTATION_FREQUENCY>\n";
	    print OUT "<ANNOTATION>$annotation\n";
	    if ($gene_annotations ne "NONE") {
		print OUT "$gene_annotations\n";
	    }
	    print OUT "</ANNOTATION>\n";
	    print OUT "</GENE_ALTERATION>\n";
	}
    }
    print OUT "</OTHER_GENOMIC_ALTERATIONS_FOUND>\n";

    if ($NEED_JOURNAL_REFERENCES == $TRUE) {
	print OUT "<REFERENCES>\n";
	my $journal = "";
	open(J, "<$journal_references") || die "Journal references file not found: $journal_references\n";
	while(my $row = <J>) {
	    chomp $row;
	    if ($row && $row ne "") {
		$journal .= "<REFERENCE>$row</REFERENCE>\n";
	    }
	}
	close(J);
	if ($journal eq "") {
	    print OUT "<REFERENCE>   </REFERENCE>\n";
	} else {
	    print OUT "$journal\n";
	}
	print OUT "</REFERENCES>\n";
    }#end if ($NEED_JOURNAL_REFERENCES == $TRUE)


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

#    foreach my $m (@data) {
#	my $test_mut = &flattenMutation($m);
#	$test_mut .= join("\t", $m->is_vcs, $m->is_validated);
#	print "LOOK: $test_mut\n";
#    }



    my @results = ();
    if (scalar(@data) == 0) {
	my $datum = join("\t", "NONE FOUND", "NA", "NA");	
	push @results, $datum;
	return (0, \@results);
    }

    foreach my $mut (@data) {
	if ($mut->is_validated eq $is_indp_validated) {
	    my $gene = $mut->gene_name;
	    my $alteration = $mut->alteration;
	    my $gene_alteration = "$gene $alteration";
	    my $mut_freq = $mut->mutation_frequency;
	    my $annotation = "";
	    my $is_vcs = $mut->is_vcs;
#	    print ">>>>>>>>> $gene_alteration\n";
#	    print ">>>>>>>>> $is_vcs\n";
#	    if ($is_vcs eq "VCS") {
	    if ($is_vcs =~ /^VCS/) {
		# We get interpretation from vcs database
		print "IS A VCS <<<<< \n";
		my $chr = $mut->chr;
		my $ref = $mut->ref;
		my $pos = $mut->pos;
		my $alt = $mut->alt;
		print "QUERYING INTERPRETATION: $disease, $chr, $ref, $pos, $alt\n";
		$stm->execute($disease, $chr, $ref, $pos, $alt); #DISEASE=? AND CHR=? AND REF_ALLE=? AND GENOMIC_COORD=? AND ALT_ALLELE=?";
		while(my @row_array = $stm->fetchrow_array()) {
		    $annotation = $row_array[0];
		}
		if (!$annotation || $annotation eq "") {
		    $annotation = &getGenericInterpretation($gene, $generic_database, $generic_table);
		}
	    } else {
		# We get from generic database
		$annotation = &getGenericInterpretation($gene, $generic_database, $generic_table);
	    }


	    my $for_txt_file = join(":", $gene, $alteration, $mut_freq);

	    my @gene_annotations = ();
	    if ($mut->cosmic_id && $mut->cosmic_id =~ /COSM/) {
		push @gene_annotations, "COSMIC: " . $mut->cosmic_id;
		$for_txt_file .= ":" . $mut->cosmic_id;
	    } else { $for_txt_file .= ":NA"; }

	    if ($mut->dbsnp && $mut->dbsnp =~ /rs/) {
		push @gene_annotations, "dbSNP: " . $mut->dbsnp;
		$for_txt_file .= ":" . $mut->dbsnp;
	    } else { $for_txt_file .= ":NA"; }
	    if ($mut->esp && $mut->esp =~ /ESP/) {
		push @gene_annotations, "ESP: " . $mut->esp;
		$for_txt_file .= ":" . $mut->esp;
	    } else { $for_txt_file .= ":NA"; }
	    if ($mut->polyphen && $mut->polyphen ne "." && $mut->polyphen ne "") {
		push @gene_annotations, "PolyPhen: " . $mut->polyphen;
		$for_txt_file .= ":" . $mut->polyphen;
	    } else { $for_txt_file .= ":NA"; }

	    my $gene_anno_line = "";
	    if (scalar(@gene_annotations) == 0) {
		$gene_anno_line = "NONE";
	    } else {
		$gene_anno_line = join(", ", @gene_annotations);
	    }

	    my $datum = join("\t", $gene_alteration, $mut_freq, $annotation, $gene_anno_line, $for_txt_file);

	    print "DATUM: $datum\n";
	    push @results, $datum;
	}
    }
    my $number_of_variants_to_report = scalar(@results);
    print "$is_indp_validated: $number_of_variants_to_report\n";
    return ($number_of_variants_to_report, \@results);
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
    if ($interpretation eq "No annotations available.") { print "GENE WITHOUT GENERIC ANNOTATIONS: $gene\n"; }
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

    print "from Step7.Generate.Clinical.Report.pl\n";
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





