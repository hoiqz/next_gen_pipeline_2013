#!/usr/bin/perl
package SNPCallingManager;
use strict;
use File::Basename;

sub new{
    my $class=shift;
    my $self={@_};
    bless $self,$class;
    return $self;
}

sub timenow{
    my $time=`date`;
    $time=~s/\n$//;
    return $time;
}

sub BEDValidation{
    my $self=shift;
    my $logger=$self->{logger};
    my $bedfile=$self->{bedfile};
    my $runtime=$self->timenow();
    my $intervalsfile=$self->{intervalsfile};

    print LOGGER "$runtime INFO: Validating BED file format\n"; 

    open(BED,"$bedfile");
    open(LOGGER,">>$logger");
    open INTERVAL,">$intervalsfile";
    while(<BED>)
    {
	next if (/^#/ || /^track/);
	chomp;
	if (/\t/)
	{
#	    print LOGGER "INFO: BEDFILE IS TAB DELIMITED FORMAT\n";
	    if($_!~/^chr/)
	    {
		next;
	    }
	    chomp;
	    my @line=split(/\t/,$_);
	    my $chrom=$line[0];
	    my $start=$line[1];
	    my $end=$line[2];
	    print INTERVAL "$chrom:$start-$end\n";
	}elsif(/:/ && /-/){
	    print LOGGER "INFO: BEDFILE IS already in CHR:start-end FORMAT\n";
	    print LOGGER "INFO: creating softlink to $intervalsfile\n";
	    system "ln -s $bedfile $intervalsfile";
	    last;
	}else{
	    print LOGGER "INFO: BEDFILE IS OF UNKNOWN FORMAT. PLEASE CHECK\n";
	    return 0;
	}
    }
    close INTERVAL;
    close LOGGER;
    return 1;
}

sub writeConfig{
    my $self=shift;
    my $platform=shift;
    my $samplename=$self->{samplename};
    my $assembly=$self->{assembly};
    my $config=$self->{config};
    my $logger=$self->{logger};
    open(LOGGER,">>$logger");
    open(SNPCONFIG,">$config");
    print LOGGER "INFO: CREATING CONFIG FILE FOR SNV CALLING AT $config\n";
    print SNPCONFIG "Library_ID\t$samplename\n";
    print SNPCONFIG "Assembly\t$assembly\n";
    print SNPCONFIG "Lib_dir\t$self->{root}\n";
    print SNPCONFIG "Output_dir\t$self->{root}\n";
    print SNPCONFIG "Log_file\t$self->{snp_calling_script_log}\n";
    print SNPCONFIG "Default_platform\t$platform\n";
    close SNPCONFIG;
    close LOGGER;
}

sub splitSNVAndIndel{
    my $self=shift;
    my $prefix=shift;
    my $snp_pileup=shift;
    my $indel_pileup=shift;
    my $splitscript=$self->{splitscript};
    my $perl=$self->{perl};
    my $logger=$self->{logger};
    open(LOGGER,">>$logger");
    my $snp_calling_pileup=$self->{snp_calling_pileup};
    my $snp_calling_root=$self->{root};
    my $waitinseconds=7200;
    print LOGGER "INFO: SPLITTING PILEUP\n";
    print LOGGER "DEBUG:$perl $splitscript $snp_calling_pileup $snp_calling_root/$prefix\n";
    my $qsubfile="$snp_calling_root/qsub_split_command.sh";
    my $splitsuccess="$snp_calling_root/SNVCalling.split.successful";
    open (COMMAND,">$qsubfile");
    print COMMAND "#!/bin/sh\n";
    print COMMAND "#\$ -S /bin/sh\n";
    print COMMAND "$perl $splitscript $snp_calling_pileup $snp_calling_root/$prefix\n";
    print COMMAND "touch $splitsuccess";
    close COMMAND;
    chmod 0775, $qsubfile;

    #system("$perl $splitscript $snp_calling_pileup $snp_calling_root/$prefix");

    #submit job to queue
    my $splitjobid=`qsub -cwd -terse $qsubfile`;

    # trace the job
    my $trace_SNV_split_success=$self->trace_completed($splitsuccess,$waitinseconds);
    if ($trace_SNV_split_success)
    {

	if((! -e $snp_pileup) && ( ! -e $indel_pileup))
	{
	    print LOGGER "ERROR: Error splitting pileup $snp_calling_pileup into $snp_pileup and $indel_pileup\n";
	    return 0;
	}else{
	    print LOGGER "INFO : Splitting pileup $snp_calling_pileup into $snp_pileup and $indel_pileup completed\n";
	    return 1;
	}
    }
    close LOGGER;
}

sub trace_completed{
    my $self=shift;
    my $trace_file=shift;
    my $wait_counter=shift;# how many seconds to wait
    my $notdone=1;
    while($notdone)
    {
	if ($wait_counter == 0)
	{
	    return 0;
	}
	if( -e $trace_file)
	{
	    # tracefile found! run completed
	    $notdone=0;
	    return 1;
	}
	else{
	    sleep(60);
	}
	$wait_counter--;
    }
}

sub runSNVCalling{
    my $self=shift;
    my $logger=$self->{logger};
    my $snp_calling_script=$self->{snp_calling_script};
    my $snp_calling_log=$self->{snp_calling_script_log};
    my $snp_calling_pileup=$self->{snp_calling_pileup};
    my $snp_calling_root=$self->{root};
    my $samplename=$self->{samplename};
    my $assembly=$self->{assembly};
    my $intervalsfile=$self->{intervalsfile};
    my $waitinseconds=86400;

    my $runtime=$self->timenow();

    print LOGGER "$runtime INFO: Clearing old success file in $snp_calling_root\n";
    my @oldsuccessfiles=<$snp_calling_root/SNVpipeline.successful>;
    foreach (@oldsuccessfiles)
    {
	unlink;
    }

    open(LOGGER,">>$logger");
    print LOGGER "$runtime INFO: SUBMITTING SNV CALLING\n";
    print "$runtime INFO: SUBMITTING SNV CALLING\n";
    print LOGGER "$runtime INFO: qsub -cwd -terse -q p_variantcalling.q $snp_calling_script $snp_calling_root $samplename $assembly $intervalsfile\n";

    # submit job to queue
    my $snpcallingjobid=`qsub -cwd -terse -q p_variantcalling.q $snp_calling_script $snp_calling_root $samplename $assembly $intervalsfile`;

    # trace the job
    my $SNVpipeline_success="$snp_calling_root/SNVpipeline.successful";
    my $trace_SNV_pipeline_success=$self->trace_completed($SNVpipeline_success,$waitinseconds);
    if ($trace_SNV_pipeline_success)
    {
	print "DEBUG: SNVpipeline.successful detected at $snp_calling_root/SNVpipeline.successful\n";
	#check the snp calling log file if it reports success
	my $pileupfound=0;
	while(! $pileupfound)
	{
	    unless (-e $snp_calling_log){
		print LOGGER "ERROR: CANT FIND LOG FILE $snp_calling_log\n";
		print "ERROR: CANT FIND LOG FILE $snp_calling_log\n";
		return 0;
	    }
	    my $checklog=`tail -1 $snp_calling_log`;
	    if(-e $snp_calling_pileup){
		if ($checklog=~/SNV pipeline completed successfully/)
		{
		    print LOGGER "INFO: SNV PIPELINE COMPLETED SUCCESSFULLY\n";
		    $pileupfound=1;
		    return 1;
		}
		elsif($checklog=~/ERROR/){
		    print LOGGER "ERROR: SNV PIPEPLINE FAILED. CHECK THE LOG AT $snp_calling_log\n";
		    return 0;
		}
		else{
		    print LOGGER "ERROR: UNCAUGHT LOG ERROR. CHECK LOG AT $snp_calling_log\n";
		    print "ERROR: UNCAUGHT LOG ERROR. CHECK LOG AT $snp_calling_log\n";
		    return 0;
		}
	    }
	    else{
		print LOGGER "ERROR: PILEUP FILE NOT FOUND. CHECK LOG AT $snp_calling_log\n";
		print "ERROR: PILEUP FILE NOT FOUND. CHECK LOG AT $snp_calling_log\n";
		return 0;
	    }
	}
    }
    else{
	print LOGGER "ERROR RUNNING SNVPIPELINE. CHECK LOG AT $snp_calling_log\n";
	print "ERROR RUNNING SNVPIPELINE. CHECK LOG AT $snp_calling_log\n";
	return 0;
    }
    close LOGGER;
}

1;
