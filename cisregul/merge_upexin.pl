#!/usr/bin/perl
# merge_upexin.pl
# given 3 separate files with upstream regions, first exons, first introns merge all these into a single upstream file for a group of genes
#
# 11-10-2003: output a list of the upstream, exon and intron length for each sequence
# 11-10-2003: make separator an optional input
# 10-01-2003: three separate information files, one for upstream, one for exon1, one for intron1. three separate indices as well
# 10-01-2003: also, separate max subseq length values for upstream, exon1 and intron1.

$CODE_DIR=$ENV{'CODEDIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";

$n_args=$#ARGV+1;
if ($n_args < 5) {
    print "usage:\n";
    print "merge_upexin.pl <ids_filename> <info_upstream_filename> <info_exon1_filename> <info_intron1_filename> <up_seq_filename> <exon1_seq_filename> <intron1_seq_filename> (<max_upstream_length> <max_exon1_length> <max_intron1_length> <print_notfound> <load2memory> <separator>)\n";
    print "\nwhere\t<ids_filename> = file with ids to search for\n";
    print "\tnote: the program searches in all columns. therefore, one should format this ids file to contain only columns with relevant information!\n";
    print "\t<info_upstream_filename> = information file for the upstream sequences (used to retrieve the line number\n";
    print "\t<info_exon1_filename> = information file for the exon1 sequences (used to retrieve the line number\n";
    print "\t<info_intron1_filename> = information file for the intron1 sequences (used to retrieve the line number\n";
    print "\t<up_seq_filename> = upstream sequences\n";
    print "\t<exon1_seq_filename> = exon 1 sequences\n";
    print "\t<intron1_seq_filename> = intron 1 sequences\n";
    print "\t<max_upstream_length> = maximum subsequence length for upstream sequences [default=5000]\n";
    print "\t<max_exon1_length> = maximum subsequence length for exon1 sequences [default=5000]\n";
    print "\t<max_intron1_length> = maximum subsequence length for intron1 sequences [default=5000]\n";
    print "\t<print_notfound> = 1 to also print the information for the entries that are not found without the sequence and 0 to skip these entries [default=1]\n";
    print "\t<load2memory> = 1 to load all the sequences to memory (speeding up the program a little bit [default=0]\n";
    print "\t<separator> = nothing if separator is 'none', otherwise use this input to separate upstream, exon and intron when they are truncated [default='nnnnnnn']\n";
    print "\n\noutput is stored in log.txt (log information) and out.txt (sequences in fasta format)\n";
    gk();
    exit;
}

$ids_filename=$ARGV[0];
$info_upstream_filename=$ARGV[1];
$info_exon1_filename=$ARGV[2];
$info_intron1_filename=$ARGV[3];
$up_seq_filename=$ARGV[4];
$exon1_seq_filename=$ARGV[5];
$intron1_seq_filename=$ARGV[6];
if ($n_args>7) {
    $max_upstream_length=$ARGV[7];
} else {
    $max_upsteram_length=5000;
}
if ($n_args>8) {
    $max_exon1_length=$ARGV[8];
} else {
    $max_exon1_length=5000;
}
if ($n_args>9) {
    $max_intron1_length=$ARGV[9];
} else {
    $max_intron1_length=5000;
}
if ($n_args>10) {
    $print_notfound=$ARGV[10];
} else {
    $print_notfound=1;
}
if ($n_args>11) {
    $load2memory=$ARGV[11];
} else {
    $load2memory=0;
}
if ($n_args>12) {
    $separator=$ARGV[12];
    if ($separator eq "none") {
	$separator="";
    }
} else {
    $separator="nnnnnnn";
}

############
# log file #
############
open (log_file,">log.txt") || die "could not open log.txt for output";
$t=get_time_string();
print log_file "merge_upexin.pl\n";
print log_file "$t\n";
print log_file "merge_upexein.pl ";
for ($i=0;$i<$n_args;$i++) {
    print log_file "$ARGV[$i] ";
}
print log_file "\n";
print log_file "ids_filename=$ids_filename\n";
print log_file "info_upstream_filename=$info_upstream_filename\n";
print log_file "info_exon1_filename=$info_exon1_filename\n";
print log_file "info_intron1_filename=$info_intron1_filename\n";
print log_file "up_seq_filename=$up_seq_filename\n";
print log_file "exon1_seq_filename=$exon1_seq_filename\n";
print log_file "intron1_seq_filename=$intron1_seq_filename\n";
print log_file "max_upstream_length=$max_upstream_length\n";
print log_file "max_exon1_length=$max_exon1_length\n";
print log_file "max_intron1_length=$max_intron1_length\n";
print log_file "print_notfound=$print_notfound\n";

###############
# get indices #
###############

# retrieve the information files
@info_upstream=read_whole_file_v4($info_upstream_filename,1,1);
@info_exon1=read_whole_file_v4($info_exon1_filename,1,1);
@info_intron1=read_whole_file_v4($info_intron1_filename,1,1);

open (ids_file,$ids_filename) || die "could not open $ids_filename";

$j=0;
while ($record=<ids_file>) {
    chomp $record;
    @splitrecord=split /\t/,$record;
    $n=$#splitrecord+1;

    # get line number for upstream sequences
    $i=0;
    $ln=-1;
    $gotit=0;
    while ( ($i<$n) & ($gotit==0) ) {
	$entry=$splitrecord[$i];
	if ($entry ne "-1") {
	    @searchstuff=grep /\b${entry}\b/,@info_upstream;
	    if (@searchstuff) {
		$n_matches=$#searchstuff+1;
		if ($n_matches > 1) {
		    print "WARNING! entry=$entry\tn_matches=$n_matches\n";
		    print log_file "WARNING! entry=$entry\tn_matches=$n_matches in upstream\tchoosing first entry\n";
		}
		$ln=get_line_number($searchstuff[0]);
		$gotit=1;
	    } else {
		print log_file "${entry}\tnot found in info_upstream\n";
		print "${entry}\tnot found in info_upstream\n";
	    }
	}
	$i++;
    }    
    $j++;                           # one more input query
    chomp $ln;
    $line_upstream[$j]=$ln;
     print "record=$record\tln=$ln...\n";

    # get line number for exon1 sequences
    $i=0;
    $ln=-1;
    $gotit=0;
    while ( ($i<$n) & ($gotit==0) ) {
	$entry=$splitrecord[$i];
	if ($entry ne "-1") {
	    @searchstuff=grep /\b${entry}\b/,@info_exon1;
	    if (@searchstuff) {
		$n_matches=$#searchstuff+1;
		if ($n_matches > 1) {
		    print "WARNING! entry=$entry\tn_matches=$n_matches\n";
		    print log_file "WARNING! entry=$entry\tn_matches=$n_matches in exon1\tchoosing first entry\n";
		}
		$ln=get_line_number($searchstuff[0]);
		$gotit=1;
	    } else {
		print log_file "${entry}\tnot found in info_exon1\n";
		print "${entry}\tnot found in info_exon1\n";
	    }
	}
	$i++;
    }    
    chomp $ln;
    $line_exon1[$j]=$ln;

    # get line number for intron1 sequences
    $i=0;
    $ln=-1;
    $gotit=0;
    while ( ($i<$n) & ($gotit==0) ) {
	$entry=$splitrecord[$i];
	if ($entry ne "-1") {
	    @searchstuff=grep /\b${entry}\b/,@info_intron1;
	    if (@searchstuff) {
		$n_matches=$#searchstuff+1;
		if ($n_matches > 1) {
		    print "WARNING! entry=$entry\tn_matches=$n_matches\n";
		    print log_file "WARNING! entry=$entry\tn_matches=$n_matches in intron1\tchoosing first entry\n";
		}
		$ln=get_line_number($searchstuff[0]);
		$gotit=1;
	    } else {
		print log_file "${entry}\tnot found in info_intron1\n";
		print "${entry}\tnot found in info_intron1\n";
	    }
	}
	$i++;
    }    
    chomp $ln;
    $line_intron1[$j]=$ln;
}
close (ids_file);

print log_file "read $j information lines\n";
print log_file "indices:\n";
print log_file "upstream\texon1\tintron1\n";
for ($i=1;$i<=$j;$i++) {
    print log_file "$line_upstream[$i]\t$line_exon1[$i]\t$line_intron1[$i]\n";
    print "$line_upstream[$i]\t$line_exon1[$i]\t$line_intron1[$i]\n";
}


if ($load2memory eq 1) {
    #####################
    # reading sequences #
    #####################

    print "reading upstream sequences...\n";   
    @upstream_sequences=read_whole_fasta($up_seq_filename);
    print "reading exon1 sequences...\n";
    @exon1_sequences=read_whole_fasta($exon1_seq_filename);
    print "reading intron1 sequences...\n";
    @intron1_sequences=read_whole_fasta($intron1_seq_filename);
}

$n_hits=0;
$n_miss=0;

print "retrieving all sequences...\n";
open (output_file,">merge_upexin.out.txt") || die "could not open merge_upexin.out.txt for output";
open (lengths_file,">merge_upexin.lengths.txt") || die "could not open merge_upexin.lengths.txt for output";

for ($i=1;$i<=$j;$i++) {

    $curr_upstream_index=$line_upstream[$i];
    print "i=$i\tcurr_upstream_index=$curr_upstream_index\n";
    if ($curr_upstream_index > -1) {
	if ($load2memory eq 1) {
	    $info_upstream=$upstream_sequences[$curr_upstream_index][1];
	    $seq_upstream=$upstream_sequences[$curr_upstream_index][2];
	} else {
	    ($seq_upstream,$info_upstream,$errorcode)=extract_seq_from_fasta($up_seq_filename,$curr_upstream_index,0);
	}
	$orig_ul=length($seq_upstream);           # original upstream sequence length
	if ($orig_ul>$max_upstream_length) {
	    $seq=extract_substr_threeprime($seq_upstream,$max_upstream_length);
	    $seq_upstream="${separator}${seq}";
	    print log_file "i=$i\tupstream\torig_ul=$orig_ul\tshortened to $max_upstream_length\n";
	    $ul=$max_upstream_length;
	} else {
	    $ul=$orig_ul;
	}
    } else {
	$info_upstream="none";
	$seq_upstream="";
	print log_file "$i\tcurr_upstream_index=$curr_upstream_index\n";
	$ul=0;
    }

    $curr_exon1_index=$line_exon1[$i];
    print "i=$i\tcurr_exon1_index=$curr_exon1_index\n";
    if ($curr_exon1_index > -1) {
	if ($load2memory eq 1) {
	    $info_exon1=$exon1_sequences[$curr_exon1_index][1];
	    $seq_exon1=$exon1_sequences[$curr_exon1_index][2];
	} else {
	    ($seq_exon1,$info_exon1,$errorcode)=extract_seq_from_fasta($exon1_seq_filename,$curr_exon1_index,0);
	}
	$orig_el=length($seq_exon1);           # original exon1 sequence length
	if ($orig_el>$max_exon1_length) {
	    $seq=substr($seq_exon1,0,$max_exon1_length);          # extract from the 5' end
	    $seq_exon1="${seq}${separator}";
	    print log_file "i=$i\texon1\torig_el=$orig_el\tshortened to $max_exon1_length\n";
	    $el=$max_exon1_length;
	} else {
	    $el=$orig_el;
	}
    } else {
	$info_exon1="none";
	$seq_exon1="";
	print log_file "$i\tcurr_exon1_index=$curr_exon1_index\n";
	$el=0;
    }

    $curr_intron1_index=$line_intron1[$i];
    print "i=$i\tcurr_intron1_index=$curr_intron1_index\n";
    if ($curr_intron1_index > -1) {
	if ($load2memory eq 1) {
	    $info_intron1=$intron1_sequences[$curr_intron1_index][1];
	    $seq_intron1=$intron1_sequences[$curr_intron1_index][2];
	} else {
	    ($seq_intron1,$info_intron1,$errorcode)=extract_seq_from_fasta($intron1_seq_filename,$curr_intron1_index,0);
	}
	$orig_il=length($seq_intron1);           # original intron1 sequence length
	if ($orig_il>$max_intron1_length) {
	    $seq=substr($seq_intron1,0,$max_intron1_length);
	    $seq_intron1="${seq}${separator}";
	    print log_file "i=$i\tintron1\torig_il=$orig_il\tshortened to $max_intron1_length\n";
	    $il=$max_intron1_length;
	} else {
	    $il=$orig_il;
	}
    } else {
	$info_intron1="none";
	$seq_intron1="";
	print log_file "$i\tcurr_intron1_index=$curr_intron1_index\n";
	$il=0;
    }

    if ( ($ul>0) | ($el>0) | ($il>0) ) {
	$s="${seq_upstream}${seq_exon1}${seq_intron1}";
	$seq=singleline2multipleline($s,80);
	@split_information=split /\t/,$info_upstream;
	print output_file ">$i\t$split_information[1]\t$split_information[2]\t$split_information[3]\t$split_information[4]\t$ul\t$orig_ul\t$el\t$orig_el\t$il\t$orig_il\n";
	print output_file "${seq}\n";
	print lengths_file "${i}\t${ul}\t${el}\t${il}\n";
	$n_hits++;
    } else {
	if ($print_notfound eq 1) {
	    print output_file ">$i\t$curr_index_upstream\t$curr_index_exon1\t$curr_index_intron1\n";
	}
	$n_miss++;
    }
}
close (output_file);
close (lengths_file);
print "processed $j entries\n";
print "n_hits=$n_hits\n";
print "n_miss=$n_miss\n";
print log_file "processed $j entries\n";
print log_file "n_hits=$n_hits\n";
print log_file "n_miss=$n_miss\n";
close (log_file);

