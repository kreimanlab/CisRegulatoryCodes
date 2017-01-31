#!/usr/bin/perl
# retrieve_upstream_v3.pl
# retrieve upstream, exon and intron
# this is a rather different version from retrieve_upstream.pl (actually, it comes from merge_upexin.pl; in contrast to merge_upexin.pl here the upstream, exon and intron files are taken from the ncbi retrieved sequences and they are not input)

# new in version 3
# 11-17-03 added fasta to the sequence filenames
# 11-17-03 exon -> exon1 and intron -> intron1 in sequence file names
# output a list of the upstream, exon and intron length for each sequence
# for the exon and intron, shorten from the 5' end (if we need to shorten at all)
# accept a separator as input (therefore allowing no separator if separator is none)
# add system_status as output

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
$prog_version=3;
$program_name="retrieve_upstream_v${prog_version}.pl";
$system_status=0;

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "${program_name}\n";
    print "usage:\n";
    print "${program_name} <ids_filename> <species> (<max_upstream> <max_exon> <max_intron> <print_notfound> <separator>)\n";
    print "\nwhere\t<ids_filename> = file with ids to search for\n";
    print "\t<species> = hs/mm\n";
    print "\t<max_upstream/exon/intron> = maximum subsequence length for upstream, exon and intron sequences respectively [default=10000]\n";
    print "\t<print_notfound> = 1 to also print the information for the entries that are not found without the sequence and 0 to skip these entries [default=1]\n";
    print "\t<separator> = use this as a separator between upstream / exon / intron if the sequences had to be shortened\n";
    print "\noutput is stored in retrieve_upstream_v${prog_version}.log.txt, retrieve_upstream_v${prog_version}.out.txt and retrieve_upstream_v${prog_version}.lengths.txt\n";
    gk();
    exit;
}


$ids_filename=$ARGV[0];
$species=$ARGV[1];
if ($n_args>2) {
    $max_upstream=$ARGV[2];
} else {
    $max_upstream=10000;
}
if ($n_args>3) {
    $max_exon=$ARGV[3];
} else {
    $max_exon=10000;
}
if ($n_args>4) {
    $max_intron=$ARGV[4];
} else {
    $max_intron=10000;
}
if ($n_args>5) {
    $print_notfound=$ARGV[5];
} else {
    $print_notfound=1;
}
if ($n_args>6) {
    $separator=$ARGV[6];
    if ($separator eq "none") {
	$separator="";
    }
} else {
    $separator="";
}

# open log file
open (log_file,">retrieve_upstream_v${prog_version}.log.txt") || die "could not open log.txt for output";
print log_file "$program_name\n";
$t=get_time_string();
print log_file "$t\n";
print log_file "${program_name} ${ids_filename} ${species} ${max_upstream} ${max_exon} ${max_intron} ${print_notfound} ${separator}\n";

# get chromosome
$filename="${DATA_DIR}/db/upstream_sequences/${species}/lists/exon_list.txt";
print log_file "list filename=$filename\n";
@info=read_whole_file_v4($filename,1,1);
print log_file "info[1]=$info[1]\n";

print log_file "ids_filename=$ids_filename\n";
open (ids_file,$ids_filename) || die "could not open $ids_filename";

$j=0;
while ($record=<ids_file>) {
    chomp $record;
    @splitrecord=split /\t/,$record;
    $gotit=0;
    $n=$#splitrecord+1;
    $i=0;
    $curr_chr=-1;
    while ( ($i<$n) & ($gotit==0) ) {
	$entry=$splitrecord[$i];
	if ($entry ne "-1") {
	    if ($entry =~ /^ll/) {
		$entry=~s /ll\./LocusID:/;
	    }
	    if ($entry =~ /^\d/) {
		# if it starts directly with a number, assume it is a locus id
		$entry="LocusID:${entry}";
	    }
	    @searchstuff=grep /\b${entry}\b/,@info;
	    if (@searchstuff) {
		$n_matches=$#searchstuff+1;
		if ($n_matches > 1) {
		    print "WARNING! entry=$entry\tn_matches=$n_matches\n";
		    print log_file "WARNING! entry=$entry\tn_matches=$n_matches\tchoosing first entry\n";
		}
		@splitarray=split /\|/,$searchstuff[0];
		$curr_chr=$splitarray[1];
		$curr_chr=~s/\s+//g;
		$curr_chr=~s/${species}_chr//;
		$curr_chr=~s/\.gbk//;
		$gotit=1;
	    } else {
		print log_file "${entry}\tnot found in info\n";
		print "${entry}\tnot found in info\n";
	    }
	}
	$i++;
    }
    $j++;
    $chr[$j]=$curr_chr;
}
$n_entries=$j;
close (ids_file);
print log_file "read $j information lines\n";
print log_file "n_entries=$n_entries\n";
print log_file "chromosomes:\t";
for ($i=1;$i<=$n_entries;$i++) {
    print log_file "$chr[$i]\t";
}
print log_file "\n";

# get the indices in each chromosome file
open (ids_file,$ids_filename) || die "could not open $ids_filename for reading";
$j=1;
while ($record=<ids_file>) {
    chomp $record;
    $curr_exon_index=-1;
    $curr_intron_index=-1;
    $curr_upstream_index=-1;

    $curr_chr=$chr[$j];
    if ($curr_chr ne "-1") {
	@splitrecord=split /\t/,$record;
	$n=$#splitrecord+1;

	$i=0;
	$gotit=0;
	$filename="${DATA_DIR}/db/upstream_sequences/${species}/lists/${species}_chr${curr_chr}.exon_list.txt";
	print log_file "filename=$filename\n";
	@info=read_whole_file_v4($filename,1,1);
	while ( ($i<$n) & ($gotit==0) ) {
	    $entry=$splitrecord[$i];
	    if ($entry ne "-1") {
		if ($entry =~ /^ll/) {
		    $entry=~s /ll\./LocusID:/;
		}
		if ($entry =~ /^\d/) {
		    # if it starts with a number, assume it is a locus id 
		    $entry="LocusID:${entry}";
		}
		($n_matches,$curr_exon_index,$gotit)=get_index($entry,@info);
		print log_file "\t$n_matches\t$curr_exon_index\t$entry\n";
	    } else {
		print log_file "${entry}\tnot found in info\n";
		print "${entry}\tnot found in info\n";
	    }
	    $i++;
	}

	$i=0;
	$gotit=0;
	$filename="${DATA_DIR}/db/upstream_sequences/${species}/lists/${species}_chr${curr_chr}.intron_list.txt";
	print log_file "filename=$filename\n";
	@info=read_whole_file_v4($filename,1,1);

	while ( ($i<$n) & ($gotit==0) ) {
	    $entry=$splitrecord[$i];
	    if ($entry ne "-1") {
		if ($entry =~ /^ll/) {
		    $entry=~s /ll\./LocusID:/;
		}
		if ($entry =~ /^\d/) {
		    # if it starts with a number, assume it is a locus id 
		    $entry="LocusID:${entry}";
		}
		($n_matches,$curr_intron_index,$gotit)=get_index($entry,@info);
		print log_file "\t$n_matches\t$curr_intron_index\t$entry\n";
	    } else {
		print log_file "${entry}\tnot found in info\n";
		print "${entry}\tnot found in info\n";
	    }
	    $i++;
	}
	
	$i=0;
	$gotit=0;
	$curr_upstream_index=-1;
	$filename="${DATA_DIR}/db/upstream_sequences/${species}/lists/${species}_chr${curr_chr}.upstream_list.txt";
	print "filename=$filename\n";
	@info=read_whole_file_v4($filename,1,1);
	while ( ($i<$n) & ($gotit==0) ) {
	    $entry=$splitrecord[$i];
	    if ($entry ne "-1") {
		if ($entry =~ /^ll/) {
		    $entry=~s /ll\./LocusID:/;
		}
		if ($entry =~ /^\d/) {
		    # if it starts with a number, assume it is a locus id 
		    $entry="LocusID:${entry}";
		}
		($n_matches,$curr_upstream_index,$gotit)=get_index($entry,@info);
		print log_file "\t$n_matches\t$curr_upstream_index\t$entry\n";
	    } else {
		print log_file "${entry}\tnot found in info\n";
		print "${entry}\tnot found in info\n";
	    }
	    $i++;
	}	
    }
    $exon_indices[$j]=$curr_exon_index;
    $intron_indices[$j]=$curr_intron_index;
    $upstream_indices[$j]=$curr_upstream_index;
    $j++;
}
close (ids_file);
print log_file "indices:\n";
print log_file "i\tchr\tupstream_index\texon_index\tintron_index\n";
for ($i=1;$i<=$n_entries;$i++) {
    print log_file "$i\t$chr[$i]\t$upstream_indices[$i]\t$exon_indices[$i]\t$intron_indices[$i]\n";
}
print log_file "\n";

print "retrieving sequences...\n";
open (ids_file,$ids_filename) || die "could not open $ids_filename for reading";
open (output_file,">retrieve_upstream_v${prog_version}.out.txt") || die "could not open out.txt for output";
open (lengths_file,">retrieve_upstream_v${prog_version}.lengths.txt") || die "could not open lengths.txt for output";

print log_file "curr_entry\tul_curr}\tul_orig\tel_curr\tel_orig\til_curr\til_orig\n";

for ($curr_entry=1;$curr_entry<=$n_entries;$curr_entry++) {
    $record=<ids_file>;
    $curr_chr=$chr[$curr_entry];
    print "processing i=$curr_entry\n";

    if ($curr_chr ne "-1") {

	# get upstream
	$curr_upstream_index=$upstream_indices[$curr_entry];
	if ($curr_upstream_index ne "-1") {
	    $filename="${DATA_DIR}/db/upstream_sequences/${species}/${species}_chr${curr_chr}.gbk.upstream.fasta";
	    ($entry,$info_up,$errorcode)=extract_seq_from_fasta($filename,$curr_upstream_index,0);           # returns the sequence and info (not in fasta format)
	    if ($errorcode eq 0) {
		#($upstream,$shortened)=shorten_fastaseq($entry,$max_upstream,0);
		($upstream,$shortened,$upstream_length_orig,$upstream_length_curr)=shorten_fastaseq($entry,$max_upstream,0,3,$separator);
	    }
	} else {
	    $upstream="";
	}

	# get exon1
	$curr_exon_index=$exon_indices[$curr_entry];
	if ($curr_exon_index ne "-1") {
	    $filename="${DATA_DIR}/db/upstream_sequences/${species}/${species}_chr${curr_chr}.gbk.exon1.fasta";
	    ($entry,$info_ex,$errorcode)=extract_seq_from_fasta($filename,$curr_exon_index,0);
	    if ($errorcode eq 0) {
		#($exon,$shortened)=shorten_fastaseq($entry,$max_exon,0);
		($exon,$shortened,$exon_length_orig,$exon_length_curr)=shorten_fastaseq($entry,$max_exon,0,5,$separator);
	    }
	} else {
	    $exon="";
	}

	# get intron1
	$curr_intron_index=$intron_indices[$curr_entry];
	if ($curr_intron_index ne "-1") {
	    $filename="${DATA_DIR}/db/upstream_sequences/${species}/${species}_chr${curr_chr}.gbk.intron1.fasta";
	    ($entry,$info_in,$errorcode)=extract_seq_from_fasta($filename,$curr_intron_index,0);
	    if ($errorcode eq 0) {
		#($intron,$shortened)=shorten_fastaseq($entry,$max_intron,0);
		($intron,$shortened,$intron_length_orig,$intron_length_curr)=shorten_fastaseq($entry,$max_intron,0,5,$separator);
	    }
	} else {
	    $intron="";
	}
    
	$seq="${upstream}${exon}${intron}";
	$fastaseq=singleline2multipleline($seq,80);
	$entry="${info_ex}\n$fastaseq\n";
	print output_file "$entry";
	print lengths_file "${curr_entry}\t${upstream_length_curr}\t${exon_length_curr}\t${intron_length_curr}\n";
	print log_file "${curr_entry}\t${upstream_length_curr}\t${upstream_length_orig}\t${exon_length_curr}\t${exon_length_orig}\t${intron_length_curr}\t${intron_length_orig}\n";
    } else {
	if ($print_notfound eq 1) {
	    chomp $record;
	    print output_file ">$record\n";
	    print lengths_file "${curr_entry}\t0\t0\t0\n";
	}
    }        # close check on curr_chr
}
close (output_file);
close (ids_file);
close (lengths_file);
close (log_file);
$curr_entry--;
print " ($curr_entry) ready\n";
print "system_status=$system_status\n";

# end main

sub get_index{
    # get_index
    # $id=get_index($entry,@info);

    my $entry;
    my @info;
    my @searchstuff;
    my $n_matches;
    my $id;
    my $n_matches=-1;
    my $gotit=0;

    ($entry,@info)=@_;

    @searchstuff=grep /\b${entry}\b/,@info;
    if (@searchstuff) {
	$n_matches=$#searchstuff+1;
	if ($n_matches > 1) {
	    print "WARNING! entry=$entry\tn_matches=$n_matches\n";
	    print log_file "WARNING! entry=$entry\tn_matches=$n_matches\tchoosing first entry\n";
	}
	$id=get_line_number($searchstuff[0]);
	$gotit=1;
	#print "entry=$entry\tsearchstuff=$searchstuff[0]\t$id\n";
	#exit;
    } else {
	$id=-1;
    }

    return ($n_matches,$id,$gotit);
}

sub extract_seq_from_fasta{
    # given a fasta-formatted file, extract a sequence corresponding to a given index
    # usage:
    # ($entry,$info,$errorcode)=extract_seq_from_fasta($filename,$id,$returnfasta);
    # if returnfasta=1 then return fasta formattted sequence, else return sequence

    my $filename=$_[0];
    my $id=$_[1];
    my $entry="";
    my $input_file;
    my $k;
    my $seq;
    my $errorcode=0;
    my $info;
    my $seq;
    my $returnfasta=$_[2];

    $input_file=open_read_file($filename);
    $k=0;
    while ($k<$id) {
	$entry=read_entry_beginning_with($input_file,">");
	if (!$entry) {
	    last;
	}
	$k++;
    }
    if ($k<$id) {
	print "error! end of file before reaching id=$id\n";
	$entry="";
	$errorcode=1;
    } 
    
    if ($return_fasta eq 1) {
	$info="";
	return ($entry,$info,$errorcode);
    } else {
	($info,$seq)=split_infoseq($entry);
	return ($seq,$info,$errorcode);
    }
}

sub shorten_fastaseq{
    # shorten a fasta sequence (with info and all)
    # usage:
    # ($newentry,$shortened,$seq_length_orig,$seq_length_curr)=shorten_fastaseq($entry,$max_seq_length,$fastaseq,$threeorfive,$separator);
    # if $fastaseq=1, receive and return fasta formatted entry
    # if $threeorfive is 3, then keep the 3' end, if it is 5, then keep the 5' end
    # return shortened=1 if the sequence was shortened and 0 otherwise
    # return also the original and current sequence lengths
    # 
    # 10-27-2003: solved a bug (the last argument was not being read!)
    # 10-27-2003: use separator as input to this program
    # 10-27-2003: add threeorfive to select which end to shorten

    my $entry=$_[0];
    my $newentry;
    #my $seq_length;
    my $max_seq_length=$_[1];
    my $seq;
    my $info;
    my $shortened=0;
    my $fastaseq=$_[2];
    my $threeorfive=$_[3];
    my $separator=$_[4];
    my $seq_length_orig;
    my $seq_length_curr;

    if ($fastaseq eq 1) {
	($info,$seq)=split_infoseq($entry);
    } else {
	$seq=$entry;
    }
    $seq_length_orig=length($seq);
    $seq_length_curr=$seq_length_orig;
    if ($seq_length_orig>$max_seq_length) {
	if ($threeorfive eq 3) {
	    $seq=extract_substr_threeprime($seq,$max_seq_length);
	    $seq="${separator}${seq}";
	} else {
	    $seq=substr($seq,0,$max_seq_length);
	    $seq="${seq}${separator}";
	}
	$seq_length_curr=$max_seq_length;
	$shortened=1;
    }
    if ($fastaseq eq 1) {
	$newentry=singleline2multipleline($seq,80);
	$newentry=">${info}\n${newentry}";
    } else {
	$newentry=$seq;
    }

    return ($newentry,$shortened,$seq_length_orig,$seq_length_curr);
}

