#!/usr/bin/perl -w

###########################################################################
#
# We match variants called against Clinically Actionable Genes list
# at this stage, we're not concerned with whether it is 
# a VCF file or a pileUp file; we just want to restrict the file
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
my $database = $references{CLINICALLY_ACTIONABLE_GENES_DB};
my $table = $references{CLINICALLY_ACTIONABLE_GENES_TABLE};

# INPUT VARIANT FILE
my $pileup_dir = $subdirs{PILEUP_DIR};
my $pileup_file = "$rootname" . "_filtered.txt";
my $inputFile = "$subrootdir/$pileup_dir/$pileup_file";

if ($userConfig{VEP_INFILE_FORMAT} =~ /vcf/i) {
    $inputFile = $userConfig{VCF_INFILE}; # VCF file provided by, say MiSeq instrument.
}
if (! -e $inputFile) {
    die "Unable to proceed, input file $inputFile does not exist.\n";
}

# CREATE OUTPUT DIRECTORY
my $outsubdir = $subdirs{LIMITED_DIR};
my $outdir = "$subrootdir/$outsubdir";
if (! -d $outdir) {
    system("mkdir -p $outdir");
}


# WE WANT 2 OUTPUT FILES
my $cagOutfile = "$outdir/$rootname.limited.to.clinically.actionable.genes";
my $nonCAGOutfile = "$outdir/$rootname.not.in.clinically.actionable.genes";

my $regions_of_interest_href = &getClinicallyActionableGeneRegions($disease, $database, $table);

&filterByClinicallyActionableGeneLoc($regions_of_interest_href, $inputFile, $cagOutfile, $nonCAGOutfile);
print "Completed: $cagOutfile\n$nonCAGOutfile\n";

sub filterByClinicallyActionableGeneLoc() {
    my ($regions_of_interest_href, $inputFile, $cagOutfile, $nonCAGOutfile) = @_;
    my %cagRegions = %{$regions_of_interest_href};
    
    open(CAGOUT, ">$cagOutfile") || die "Unable to open $cagOutfile for writing.\n";
    open(NONCAGOUT, ">$nonCAGOutfile") || die "Unable to open $nonCAGOutfile for writing.\n";
    open(IN, "<$inputFile") || die "Unable to open $inputFile\n"; 
    while(my $row = <IN>) {
	chomp $row;
	if ($row =~ /^\#/) {
	    print CAGOUT "$row\n";
	    print NONCAGOUT "$row\n";
	} else {
	    my ($chr, $pos, @others) = split(/\t/, $row);
	    my $aref = $cagRegions{$chr};
	    if (defined($aref)) { # Some chromosomes are like 'chr8_xxxx_random', also some chromosomes do NOT contain any of the disease genes
		my @regions = @{$aref};
		my $isInCAG = $FALSE;
		foreach my $reg (@regions) {
		    my ($start, $end) = split(/,/, $reg);
		    if ($start <= $pos && $pos <= $end) {
			$isInCAG = $TRUE;
			last; # short-circuit because already ascertained that it is in CAG region, so don't have to continue looping.
		    }
		}
		if ($isInCAG == $TRUE) {
		    print CAGOUT "$row\n";
		} else {
		    print NONCAGOUT "$row\n";
		}
	    }
	}
    }
    close(IN);
    close(NONCAGOUT);
    close(CAGOUT);

}#end filterByClinicallyActionableGeneLoc


sub getClinicallyActionableGeneRegions() {
    my ($disease, $database, $table) = @_;
    my %results = ();

    # Prepare database 
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$database", "", "", { RaiseError => 1, AutoCommit => 1 } );
    $dbh->do('PRAGMA synchronous=1');
    $dbh->do('PRAGMA cache_size=4000');

    my $query = "SELECT CHR, CDS_START, CDS_END FROM $table WHERE " . $disease . "=1";

    my $stm = $dbh->prepare($query);
    $stm->execute();
    while(my @row_array = $stm->fetchrow_array()) {
	my ($chr,$cds_start,$cds_end) = @row_array;
	if (!defined($results{$chr})) {
	    my @array = ();
	    $results{$chr} = \@array;
	}
	my $aref = $results{$chr};
	my @a = @{$aref};
	push @a, "$cds_start,$cds_end";
	$results{$chr} = \@a;
    }#end while
    $stm->finish();
    $dbh->disconnect();
    return \%results;
}#end getClinicallyActionableGeneRegions


__END__

