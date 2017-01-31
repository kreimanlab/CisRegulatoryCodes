#!/usr/bin/perl
# preporth4alignace_v3.pl
# prepare results of bl2seq or bayesalign (comparison of sequences between two species) for alignace, motifSampler (a motif search program)
# uses a fasta formatted sequence file and the output of bayesalign_2sets_v1.pl or bl2seq_2sets_v2.pl

# format for output of bl2seq_2sets_v2.pl
# 1            2             3       4               5                 6                7                 8                  9                   10          11
# hit_number   expectation   score   identities_nt   identities_perc   query_pos_init   query_pos_finit   subject_pos_init   subject_pos_finit   seq_query   seq_subject

# new in version 3
# 02-17-2004: option of including the sequences with no orthologous gene

# 10-22-2003: option of extracting the whole conserved segment (regardless of the appearance of "Ns")
# 04-03-2003: allow to extract sequences from the second species
# make it compatible with linux
# 06-19-2003: process new output from bl2seq_blast2sets_v2.pl and use the sequence output from the blast

$CODE_DIR=$ENV{'CODEDIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";
require "${CODE_DIR}/perl/lib/regul_search_methods_v1.pl";
require "${CODE_DIR}/perl/lib/parser_methods_v1.pl";

$prog_version=3;
$last_modified="02_17_2004";
$species=1;                             # which of the two entries in the bl2seq comparison to use
$threshold_indentities_perc=60;         # below this threshold, report the bl2seq output sequences in the orth2.txt file as well
$max_n_orth3=20;                        # in orth3.txt remove segments with more than this number of consecutive ns

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "preporth4alignace_v${prog_version}.pl\n";
    print "usage:\n";
    print "\npreporth4alignace.pl <sequence_filename> <bl2seq_results_filename> (<species> <map_filename> <include_noorth>)\n";
    print "<include_noorth>\t1 to include the genes with no orthologous entry in the second species and 0 to skip these\n";
    print "\noutput\n";
    print "\torth1.txt\tsequences using the bl2seq output sequences (this includes gaps, etc.)\n";
    print "\torth2.txt\tsequences using the whole homology segment\n";
    print "\torth3.txt\tthis is based on orth2.txt; in order to try to speed up the motif search program, remove long sequences with ns\n";
    print "\tpreporth4alignace.log.txt\tlog file\n";
    print "\nthe position indices in columns 6 and 7 are used for species = 1 [default] and positions 8 and 9 are used for species = 2\n";
    print "if species=2 we need the map_filename mapping the transcripts in species 1 to transcripts in species 2 (if there is no input map, the program assumes that transcript i in species 1 maps onto transcript i in species 2)\n";
    gk();
    exit;
}

$seq_filename=$ARGV[0];
$bl2seq_results_filename=$ARGV[1];
if ($n_args >= 3) {
    $species=$ARGV[2];
}

open (log_file,">preporth4alignace.log.txt") || die "could not open preporth4alignace.log.txt for writing";
print log_file "preporth4alignace_v${prog_version}.pl\n";
$t=get_time_string();
print log_file "last modified=${last_modified}\n";
print log_file "$t\n";
print log_file "seq_filename=$seq_filename\n";
print log_file "bl2seq_results_filename=$bl2seq_results_filename\n";
print log_file "species=$species\n";
$n_transcripts=fasta_line_count($seq_filename);
print log_file "n_transcritps=$n_transcripts\n";
if ($species eq 1) {
    $lim1=5;          # column for the start position of the similar sequence in the first species
    $lim2=6;          # "              end                                                  
    $seq_col=9;       # column for the sequence in the first species
} else {
    $lim1=7;          # column for the start position of the similar sequence in the second species
    $lim2=8;          # "              end
    $seq_col=10;      # column for the sequence in the second species
    
    if ($n_args >= 4) {
	$map_filename=$ARGV[3];
    } else {
	$map_filename="nomap";
    }
    if ($map_filename ne "nomap") {
	print log_file "reading map from $map_filename\n";
	open (map_file,$map_filename) || die "could not open $map_filename for reading";
	$i=0;
	while ($record = <map_file>) {
	    chomp $record;
	    @splitarray=split /\t/,$record;
	    $i++;
	    $map[$i][1]=$splitarray[0];
	    $map[$i][2]=$splitarray[1];
	}
	close (map_file);
    } else {
	for ($i=1;$i<=$n_transcripts;$i++) {
	    $map[$i][1]=$i;
	    $map[$i][2]=$i;
	}
    }
}
if ($n_args>=5) {
    $include_noorth=$ARGV[4];
} else {
    $include_noorth=0;
}
print log_file "include_noorth=${include_noorth}\n";

# read sequences
print log_file "reading sequences from $seq_filename\n";
$seq[0]=$seq_filename;
$info[0]=$seq_filename;
$seq_lengths[0]=$seq_filename;
$seq_file=open_read_file($seq_filename);
$i=0;
while ($i<1000000) {
    $i++;
    $entry=read_entry_beginning_with($seq_file,">");
    if (!$entry) {
	print "end of $seq_filename file!\n";
	last;
    }
    ($curr_info,$curr_seq)=split_infoseq($entry);
    $info[$i]=$curr_info;
    $seq_length=length($curr_seq);
    $seq_lengths[$i]=$seq_length;
    $seq[$i]=$curr_seq;
}
close($seq_file);
print log_file "\tread $i sequences\n";

# prepare the segment of ns for orth3
$orth3_segment="";
for ($i=1;$i<=$max_n_orth3;$i++) {                        
    $orth3_segment.="n";
}

$output_filename1=">orth1.txt";
$output_filename2=">orth2.txt";
$output_filename3=">orth3.txt";
open (output_file1,$output_filename1) || die "could not open $output_filename1";
open (output_file2,$output_filename2) || die "could not open $output_filename2";
open (output_file3,$output_filename3) || die "could not open $output_filename3";

for ($i=1;$i<=$n_transcritps;$i++) {
    $n1[$i][0]=0;$n1[$i][1]=0;$n1[$i][2]=0;$n1[$i][3]=0;
    $n2[$i][0]=0;$n2[$i][1]=0;$n2[$i][2]=0;$n2[$i][3]=0;
    $n3[$i][0]=0;$n3[$i][1]=0;$n3[$i][2]=0;$n3[$i][3]=0;
}

$n_nucleotides_all=0;
$bl2seq_file=open_read_file($bl2seq_results_filename);
$i=0;
$processed=0;
$prev_seq_index=0;
while ($i<1000000) {
#while ($i<20) {
    $entry=read_entry_beginning_with($bl2seq_file,">");
    if (!$entry) {
	print "end of $bl2seq_results_filename file\n";
	last;
    }
    if ($entry =~ /Locus/) {
	$processed++;
    }
    @entry_array=split /\n/,$entry;
    $n_entry_array=$#entry_array+1;
    if ($n_entry_array>1) {
	#$i++;
	print log_file "processing i=$i\tn_entry_array=$n_entry_array\n";
	$seq_index=get_line_number($entry_array[0]);                            # this corresponds to the index of the sequence in species 1
	$seq_index=~s/\>//;
	if ( ($seq_index-$prev_seq_index)>1)  {
	    print "seq_index=$seq_index\tprev_seq_index=$prev_seq_index\n";
	    # some entry was skipped (either because no orthologous gene was available or there was no match in the sequence)
	    if ($include_noorth eq 1) {
		for ($j=($prev_seq_index+1);$j<$seq_index;$j++) {
		    $i++;
		    $curr_seq_index=$j;
		    print "\tadding curr_seq_index=$curr_seq_index\n";
		    print log_file "$i sequence index\t>${curr_seq_index}\n";
		    if ($species eq 2) {
			$seq_index2=$map[$curr_seq_index][2];
			$curr_seq_index=$seq_index2;
		    }
		    $curr_seq=$seq[$curr_seq_index];                                        # current sequence (full)
		    $curr_seq_length=length($curr_seq);                                     # length of current sequence
		    $n_nucleotides_all+=$curr_seq_length;                                   # add to the total number of nucleotides

		    print output_file1 "$info[$curr_seq_index]\n";
		    print output_file2 "$info[$curr_seq_index]\n";
		    print output_file3 "$info[$curr_seq_index]\n";
		    $n1[$j][0]=$curr_seq_length;$n1[$j][1]=length($curr_seq);
		    $n2[$j][0]=$curr_seq_length;$n2[$j][1]=length($curr_seq);
		    $n3[$j][0]=$curr_seq_length;$n3[$j][1]=length($curr_seq);
		    $temp_seq1=$curr_seq;$temp_seq1=~s/n//g;
		    $n1[$j][2]=length($temp_seq1);$n1[$j][3]=$n1[$j][1]-$n1[$j][2];
		    $n2[$j][2]=length($temp_seq1);$n2[$j][3]=$n2[$j][1]-$n2[$j][2];
		    $n3[$j][2]=length($temp_seq1);$n3[$j][3]=$n3[$j][1]-$n3[$j][2];
		    $new_seq1=singleline2multipleline($curr_seq,80);
		    $new_seq1=lc($new_seq1);
		    print output_file1 "$new_seq1\n";
		    print output_file2 "$new_seq1\n";
		    print output_file3 "$new_seq1\n";
		}
	    }
	}
	$prev_seq_index=$seq_index;
	$i++;
	print log_file "$i sequence index\t>${seq_index}\n";

	if ($species eq 2) {
	    $seq_index2=$map[$seq_index][2];
	    $seq_index=$seq_index2;
	}
	$curr_seq=$seq[$seq_index];                                             # current sequence (full)
	$curr_seq_length=length($curr_seq);                                     # length of current sequence
	$n_nucleotides_all+=$curr_seq_length;                                   # add to the total number of nucleotides

	# initialize temp_seq to all "n"
	$temp_seq1="";
	$temp_seq2="";
	$seq_length=$seq_lengths[$seq_index];
	print log_file "\tseq_length=$seq_length\n";
	for ($j=1;$j<=$seq_length;$j++) {
	    $temp_seq1.="n";
	    $temp_seq2.="n";
	}

	for ($j=1;$j<$n_entry_array;$j++) {
	    $new_seq1="";
	    $new_seq2="";
	    @splitarray=split /\t/,$entry_array[$j];
	    $p1=$splitarray[$lim1];                                             # initial position for current blast hit
	    $p2=$splitarray[$lim2];                                             # end position for current blast hit
	    if ($p1>$p2) {                                                      # minus strand hit
		$k=$p2;                                                         # $k is > 0 for the - strand
		$p2=$p1;
		$p1=$k;
	    } else {
		$k=-1;                                                          # $k is < 0 for the + strand
	    }
	    $p1--;
	    
	    # get sequence up to the match region
	    $temp_substr1=substr($temp_seq1,0,$p1);
	    $new_seq1.="$temp_substr1";
	    $temp_substr2=substr($temp_seq2,0,$p1);
	    $new_seq2.="$temp_substr2";

	    # get sequence in the match region
	    $temp_substr1=$splitarray[$seq_col];                                # get bl2seq output for sequence
	    if ($k>0) {
		$temp_substr1=revcomplement($temp_substr1);
	    }
	    $new_seq1.="$temp_substr1";

	    $identities_perc=$splitarray[4];
	    if ($identities_perc>$threshold_identities_perc) {
		$temp_substr2=substr($curr_seq,$p1,$p2-$p1);                        # get the actual sequence
	    } else {
		$temp_substr2=$temp_substr1;
	    }
	    $new_seq2.="$temp_substr2";

	    # get sequence after the match region
	    $ii=$p2;
	    $il=$seq_length-$p2;
	    $temp_substr1=substr($temp_seq1,$ii,$il);
	    if (length($temp_substr1)<0) {
		print "$i\t$j\n";
		print "$ii\t$il\n";
		exit;
	    }
	    $new_seq1.="$temp_substr1";
	    $temp_seq1=$new_seq1;

	    $temp_substr2=substr($temp_seq2,$ii,$il);
	    if (length($temp_substr2)<0) {
		print "$i\t$j\n";
		print "$ii\t$il\n";
		exit;
	    }
	    $new_seq2.="$temp_substr2";
	    $temp_seq2=$new_seq2;
	}
	print output_file1 "$info[$seq_index]\n";
	$n1[$i][0]=$seq_length;
	$n1[$i][1]=length($new_seq1);
	$temp_seq1=$new_seq1;
	$temp_seq1=~s/n//g;
	$n1[$i][2]=length($temp_seq1);
	$n1[$i][3]=$n1[$i][1]-$n1[$i][2];
	$new_seq1=singleline2multipleline($new_seq1,80);
	$new_seq1=lc($new_seq1);
	print output_file1 "$new_seq1\n";

	print output_file2 "$info[$seq_index]\n";
	$new_seq2=lc($new_seq2);
	$new_seq3=$new_seq2;
	$n2[$i][0]=$seq_length;
	$n2[$i][1]=length($new_seq2);
	$temp_seq2=$new_seq2;
	$temp_seq2=~s/n//g;
	$n2[$i][2]=length($temp_seq2);
	$n2[$i][3]=$n2[$i][1]-$n2[$i][2];
	$new_seq2=singleline2multipleline($new_seq2,80);
	print output_file2 "$new_seq2\n";

	$new_seq3=~s/${orth3_segment}//g;
	$n3[$i][0]=$seq_length;
	$n3[$i][1]=length($new_seq3);
	$temp_seq3=$new_seq3;
	$temp_seq3=~s/n//g;
	$n3[$i][2]=length($temp_seq3);
	$n3[$i][3]=$n3[$i][1]-$n3[$i][2];
	print output_file3 "$info[$seq_index]\n";
	$new_seq3=singleline2multipleline($new_seq3,80);
	print output_file3 "$new_seq3\n";
    } else {
	if ($entry =~ /Locus/) {
	    print log_file "$processed\t$n_entry_array\tno homologous sequence\n";
	}
    }        # close check on n_entry_array>1
}            # close i loop
close ($bl2seq_file);
print log_file "processed i=$i\n";
if ($include_noorth eq 1) {
    for ($j=($i+1);$j<=$n_transcripts;$j++) {
	$curr_seq_index=$j;
	print "\tadding curr_seq_index=$curr_seq_index\n";
	print log_file "sequence index\t>${curr_seq_index}\n";
	if ($species eq 2) {
	    $seq_index2=$map[$curr_seq_index][2];
	    $curr_seq_index=$seq_index2;
	}
	$curr_seq=$seq[$curr_seq_index];                                        # current sequence (full)
	$curr_seq_length=length($curr_seq);                                     # length of current sequence
	$n_nucleotides_all+=$curr_seq_length;                                   # add to the total number of nucleotides

	print output_file1 "$info[$curr_seq_index]\n";
	print output_file2 "$info[$curr_seq_index]\n";
	print output_file3 "$info[$curr_seq_index]\n";
	$n1[$j][0]=$curr_seq_length;$n1[$j][1]=length($curr_seq);
	$n2[$j][0]=$curr_seq_length;$n2[$j][1]=length($curr_seq);
	$n3[$j][0]=$curr_seq_length;$n3[$j][1]=length($curr_seq);
	$temp_seq1=$curr_seq;$temp_seq1=~s/n//g;
	$n1[$j][2]=length($temp_seq1);$n1[$j][3]=$n1[$j][1]-$n1[$j][2];
	$n2[$j][2]=length($temp_seq1);$n2[$j][3]=$n2[$j][1]-$n2[$j][2];
	$n3[$j][2]=length($temp_seq1);$n3[$j][3]=$n3[$j][1]-$n3[$j][2];
	$new_seq1=singleline2multipleline($curr_seq,80);
	$new_seq1=lc($new_seq1);
	print output_file1 "$new_seq1\n";
	print output_file2 "$new_seq1\n";
	print output_file3 "$new_seq1\n";
    }
}
close (output_file1);
close (output_file2);
close (output_file3);

print log_file "i\tn1\t\t\t\tn2\t\t\t\tn3\n";
print log_file "\t1\t2\t3\t4\t1\t2\t3\t4\t1\t2\t3\t4\n";
for ($i=1;$i<=$n_transcripts;$i++) {
    print log_file "$i\t$n1[$i][0]\t$n1[$i][1]\t$n1[$i][2]\t$n1[$i][3]\t";
    print log_file "$n2[$i][0]\t$n2[$i][1]\t$n2[$i][2]\t$n2[$i][3]\t";
    print log_file "$n3[$i][0]\t$n3[$i][1]\t$n3[$i][2]\t$n3[$i][3]\n";
}
print log_file "n_nucleotides_all=$n_nucleotides_all\n";

close (log_file);

