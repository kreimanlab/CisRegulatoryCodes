#!/usr/bin/perl
# scan_upstream_sequences_v2.pl
# run scan_motif to scan the occurrences of motif from database
#
# gabriel kreiman
# last updated: 02-04-2004
# kreiman@mit.edu

# 02-04-2004: call scan_motif_v10.exe (solved a bug in the report of strands in scan_motif_v9.c)
# 12-07-2003: report positions with respect to TSS if we get the ul,el,il lengths as input
# 12-07-2003: call scan_motif_v9.exe 
# 12-07-2003: retrieve also the statistics and the strand information
# see full history below the code

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
$TEMP_DIR=$ENV{'TEMPDIR'};

require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";
require "${CODE_DIR}/perl/lib/math_methods.pl";

$scan_version=10;                    # c code version
$program_name="scan_upstream_sequences_v2.pl";
$last_modified="02_04_2004";
$verbose=1;                         # call scan_motif with verbose output (to be able to retrieve it using scan_output here)
$rc=1;                              # scan both strands by default
$sortpos=1;                         # sort positions by default

($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
$file_t="${Hour}${Minute}${Second}";                    # to name the temporary files

$n_args=$#ARGV+1;
if ($n_args<2) {
    print "${program_name}\n";
    print "usage:\n";
    print "perl scan_upstream_sequences_v2.pl <matrix_filename> <upstream_filename> (<threshold_column> <ul> <rc> <sortpos> <ulelil_lengths_filename>)\n";
    print "\nThis program will search for the occurrences of the corresponding matrix/matrices in the upstream sequences using scan_motif_v${scan_version}.exe\n";
    print "the threshold for the score of projecting an oligonucleotide onto the weight matrix (one for each weight matrix and for each p value) is stored in the information line within the weight matrix file\n";
    print "<threshold_column> is the column where the threshold is stored [default=38]\n";
    print "<rc> if <rc> == 1, scan both strands [default=1]\n";
    print "<sortpos> if <sortpos> == 1, then sort positions\n";
    print "\nNOTE: the upstream sequences are supposed to be in lower case\n";
    print "\noutput is stored in scan_upstream_sequences.out.txt, scan_upstream_sequences.strand.txt and scan_upstream_sequences.log.txt (in call directory)\n";
    print "\nlast modified:\t${last_modified}\n\n";
    gk();
    exit;
}
$matrix_filename=$ARGV[0];
$upstream_filename=$ARGV[1];
if ($n_args > 2) {
    $threshold_column=$ARGV[2];
} else {
    $threshold_column=38;
}
if ($n_args > 3) {
    $ul=$ARGV[3];
} else {
    $ul=100000;
}
if ($n_args > 4) {
    $rc=$ARGV[4];
} else {
    $rc=1;
}
if ($n_args > 5) {
    $sortpos=$ARGV[5];
}
if ($n_args>6) {
    $ulelil_lengths_filename=$ARGV[6];
} else {
    $ulelil_lengths_filename="none";
}
$max_ul=5000;
$max_el=5000;
$max_il=5000;

$time_string=get_time_string();

$temp_filename=">${TEMP_DIR}/scan_upstream_sequences.${file_t}.tmp";
$temp_filename_nos="${TEMP_DIR}/scan_upstream_sequences.${file_t}.tmp";
$output_filename=">scan_upstream_sequences.out.txt";
$strand_filename=">scan_upstream_sequences.strand.txt";
$log_filename=">scan_upstream_sequences.log.txt";
#$thres_index+=11; # to get the appropriate column in the thresholds file

$matrix_file=open_read_file($matrix_filename);
open (output_file,$output_filename) || die "could not open $output_filename for writing";
open (strand_file,$strand_filename) || die "could not open $strand_filename for writing";
$log_file=open_read_file($log_filename);

$headertext="% program_name=$program_name (last modified = $last_modified)\n";
$headertext.="% ${time_string}\n";
$temptext="% $program_name ";
for ($i=0;$i<$n_args;$i++) {
    $temptext.="$ARGV[$i] ";
}
$headertext.="$temptext\n";
$headertext.="% matrix_filename=$matrix_filename\n";
$headertext.="% upstream_filename=$upstream_filename\n";
$headertext.="% threshold_column=${threshold_column}\n";
$headertext.="% ul=${ul}\n";
$headertext.="% rc=${rc}\n";
$headertext.="% sortpos=${sortpos}\n";
$headertext.="% ulelil_lengths_filename\t=$ulelil_lengths_filename\n";
$headertext.="% i\tid\tmotif\tmotif_length\tn_seqs\ttranscript\tpos\ttranscript\tpos\n";

print output_file "$headertext";
print strand_file "$headertext";
print log_file "$headertext";

# read ulelil_lengths_filename
$i=0;
if ($ulelil_lengths_filename ne "none") {
    open (temp_file,$ulelil_lengths_filename) || die "could not open $ulelil_lengths_filename for reading";
    while ($record=<temp_file>) {
	chomp $record;
	@splitarray=split /\t/,$record;
	$i++;
	for ($j=0;$j<=3;$j++) {
	    $ulelil[$i][$j]=$splitarray[$j];
	}
    }
    close (temp_file);
}
if ($i<=0) {
    $ulelil_lengths_filename="none";
}

$curr_motif=0;
#while ($curr_motif<3) {
while ($curr_motif<1000000) {
    $curr_motif++;
    print "processing curr_motif=$curr_motif\n";

    # get weight matrix
    $entry=read_entry_beginning_with($matrix_file,">");
    if (! $entry) {
        print "end of matrix file\n";
        last;
    }
    @splitarray=split /\n/,$entry;
    @info=split /\t/, $splitarray[0];
    $n_seqs4motif=$info[7];
    $motif_length=$info[2];

    print $log_file "weight matrix\tcurrent motif=$curr_motif\n";
    open (temp_file,$temp_filename) || die "could not open the file $temp_filename for writing";
    for ($k=1;$k<=4;$k++) {
	$record=$splitarray[$k];
	@wm_nucleotide=split /\t/,$record;
	for ($i=0;$i<$motif_length;$i++) {
	    $xt=sprintf "%.4f",$wm_nucleotide[$i];
	    print temp_file "$xt\t";
	    print $log_file "$xt\t";
	}
	print temp_file "\n";
	print $log_file "\n";
    }
    close (temp_file);

    # get threshold
    $thres=$info[$threshold_column];
    $thres-=0.0001;                        # to eliminate any possible rounding issues
    $code_exec="${CODE_DIR}/cpp/executables/scan_motif_v${scan_version}.exe ${upstream_filename} ${temp_filename_nos} ${motif_length} ${thres} ${ul} ${rc} ${sortpos} ${verbose}";
    print $log_file "code_exec=$code_exec\n";
    @scan_output=`$code_exec`;
    #print "scan_output=@scan_output\n";
    #exit;

    print "$curr_motif\t$info[1]\t$info[2]\t";
    print output_file "$curr_motif\t$info[1]\t$info[2]\t$motif_length\t$n_seqs4motif\t";
    @searchstuff=grep /^\d/,@scan_output;                   # search for the entries starting with a number (those containing transcript and positions)
    $transcripts_with_motif[$curr_motif]=$#searchstuff;     # number of transcripts with motif
    undef (@n_occurrences);
    for ($i=0;$i<=$#searchstuff;$i++) {
	$temp=$searchstuff[$i];
	chomp $temp;
	@splitarray=split /\t/,$temp;
	$n_occurrences[$i]=$#splitarray;                    # number of occurrences per transcript
	for ($j=1;$j<=$#splitarray;$j++) {
	    $p=$splitarray[$j];
	    if ($ulelil_lengths_filename ne "none") {
		($p_tss,$errorcode)=convert_positions($splitarray[0],$p,$max_ul,$max_el,$max_il,@ulelil);
	    } else {
		$p_tss=$p;
	    }
	    #print output_file "$splitarray[0]\t$splitarray[$j]\t";
	    print output_file "$splitarray[0]\t$p_tss\t";
	}
    }
    print "\n";
    print output_file "\n";
    ($mean_occurrences,$std_occurrences,$min_occurrences,$max_occurrences)=basic_stats(0,@n_occurrences);
    $motif_mean[$curr_motif]=$mean_occurrences;
    $motif_std[$curr_motif]=$std_occurrences;
    $motif_min[$curr_motif]=$min_occurrences;
    $motif_max[$curr_motif]=$max_occurrences;
    ($n_transcripts_processed,$n_matches)=extract_text_from_array("n_transcripts_processed",@scan_output);
    ($n_transcripts_hits,$n_matches)=extract_text_from_array("n_transcripts_hits",@scan_output);
    ($n_hits,$n_matches)=extract_text_from_array("n_hits",@scan_output);
    ($total_sequence_length,$n_matches)=extract_text_from_array("total_sequence_length",@scan_output);
    ($bp_per_hit,$n_matches)=extract_text_from_array("bp_per_hit",@scan_output);
    $transcripts_processed[$curr_motif]=$n_transcripts_processed;
    $transcripts_hits[$curr_motif]=$n_transcripts_hits;
    $hits[$curr_motif]=$n_hits;
    $bp[$curr_motif]=$bp_per_hit;

    # print strand information
    print strand_file "$curr_motif\t$info[1]\t$info[2]\t$motif_length\t$n_seqs4motif\t";
    open(temp_file,"scan_motif_v${scan_version}.strand.txt") || die "could not open scan_output_v${scan_version}.strand.txt for reading";
    while ($record=<temp_file>) {
	if (is_comment($record) ne "true") {
	    chomp $record;
	    @splitarray=split /\t/,$record;
	    for ($j=1;$j<=$#splitarray;$j++) {
		print strand_file "$splitarray[0]\t$splitarray[$j]\t";
	    }
	}
    }
    close (temp_file);
    print strand_file "\n";

    $status=copyfile2handle("scan_motif_v${scan_version}.log.txt",$log_file,1);
    if ($status ne 0) {
	print "warning! i found problems while copying scan_motif_v${scan_version} log file [system status = $status]\n";
    }
    #exit;
}
close ($matrix_file);
close (output_file);
close (strand_file);

# motif stats
print $log_file "\nmotif statistics\n";
print $log_file "i\ttr_with_motif\tmean_occurrences\tstd_occurrences\tmin_occurrences\tmax_occurrences\ttr_processed\ttr_hits\thits\tbp_per_hit\n";
for ($i=1;$i<$curr_motif;$i++) {
    print $log_file "$i\t$transcripts_with_motif[$i]\t$motif_mean[$i]\t$motif_std[$i]\t$motif_min[$i]\t$motif_max[$i]\t$transcripts_processed[$i]\t$transcripts_hits[$i]\t$hits[$i]\t$bp[$i]\n";
}

close ($log_file);

# remove some of the temporary files
unlink ($temp_filename_nos);
$temp_filename="scan_motif_v${scan_version}.log.txt";
unlink ($temp_filename);
$temp_filename="scan_motif_v${scan_version}.out.txt";
unlink ($temp_filename);
$temp_filename="scan_motif_v${scan_version}.strand.txt";
unlink ($temp_filename);

# end main

sub convert_positions{
    # convert positions measured with respect to the 3' end
    # to positions measured with respect to the TSS
    # (i) x=0, position = TSS
    # (ii) x<0, position = upstream nt before TSS
    # (iii) x>0, position in exon or intron
    #
    # usage:
    # ($p_tss,$errorcode)=convert_positions($transcript_index,$position,$max_ul,$max_el,$max_il,@ulelil);

    my $transcript_index;
    my $position;
    my @ulelil;
    my $ul;
    my $el;
    my $il;
    my $max_ul;
    my $max_el;
    my $max_il;
    my $p_tss;
    my $errorcode=0;

    ($transcript_index,$position,$max_ul,$max_el,$max_il,@ulelil)=@_;
    $ul=$ulelil[$transcript_index][1];
    $el=$ulelil[$transcript_index][2];
    $il=$ulelil[$transcript_index][3];
    
    if ($ul>$max_ul) {
	$ul=$max_ul;
    }
    if ($el>$max_el) {
	$el=$max_el;
    }
    if ($il>$max_il) {
	$il=$max_il;
    }

    $p_tss=$il+$el-$position;
    if (-$ul>$p_tss) {
	$errorcode=1;
    }
    #print "position=$position\n";
    #print "transcript_index=$transcript_index\n";
    #print "ul=$ul\n";
    #print "el=$el\n";
    #print "il=$il\n";
    #print "p_tss=$p_tss\n";
    #exit;
    return ($p_tss,$errorcode);
}

# weight matrix file description
# description	col	example
# id_n	1	>1
# id	2	wasserman_sites.myf
# motif length	3	14
# mean_pos	4	nan
# std_pos	5	nan
# min_pos	6	nan
# max_pos	7	nan
# n_seqs4motif	8	48
# mean_pos_p	9	nan
# std_pos_p	10	nan
# min_pos_p	11	nan
# max_pos_p	12	nan
# n_p	13	nan
# mean_pos_n	14	nan
# std_pos_n	15	nan
# min_pos_n	16	nan
# max_pos_n	17	nan
# n_n	18	nan
# consensus	19	nnvncavbtknnnn
# i_tot	20	8.1
# i_mean	21	0.6
# maxlik_seq	22	gcagcagctgcccc
# variety	23	4
# max_score	24	9.4
# min_score	25	-2
# mean_score	26	6.8
# std_score	27	2.3
# min_possible_score	28	-18
# max_possible_score	29	11
# p5	30	2.2
# p10	31	4.4
# mc	32	
# mc_thres	33	
# p0.99	34	
# p0.999	35	
# p0.9999	36	
# p0.99999	37	
# p0.999999	38	
# threshold	39	

# 08-27-2003: call scan_motif_v8.exe
# scan_motif_v8 <sequence> <weight_matrix> <motif_length> <threshold> <max_upstream_length> <rc> <sortpos> <verbose>
# 07-23-2003: use time to name the temporary files
# 07-23-2003: added statistics
# 07-17-2003: call scan_motif_v7
# 07-17-2003: retrieve thresholds from the information line in the weight matrix file
# 06-24-2003: call scan_motif_v6 (using sequences of arbitrary length: use a very large value of ul, e.g. ul=100000 so that the sequence is never truncated)
# 05-30-2003: call scan_motif_v5 (scan both + and - strands)
# 05-06-2003: call scan_motif_v4
#    usage: scan_motif_v4 <background_sequence> <weight_matrix> <motif_length> <threshold> <n_seqs4motif> <max_upstream_length> <verbose>
# 05-06-2003: subtract 0.01 to the threshold (instead of 0.1 as we did before)
# 03-03-2003: subtracted 0.1 to the threshold to correct for the 2 decimal precision (this had caused a lot of trivial problems)
# 03-24-2003: write output to scan_upstream_sequences.out.txt
