#!/usr/bin/perl
#parser_methods_v1.pl
#methods (subs) used to parse text

# list of methods
#
# $record=read_entry_beginning_with(FILEHANDLE,$delimiter);
# @query=search_pattern($entry,$search_pattern);
# ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos)=extract_blast_results_from_entry($entry,$v,$mask_nonconserved);
# ($seq)=get_sequence_from_queryline($queryline);
# ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$seq_query,$query_pos,$subject_pos)=extract_blast_results_from_array($mask_nonconserved,@blast_output);
# $pos=get_pos_from_queryline($queryline,$querysubj);
# ($info,$seq)=split_infoseq($entry);
# $txt2=singleline2multipleline($txt1,$n);
# ($score,$expect,$identities_nt,$identities_perc,$seq)=extract_blast_results_from_file($blast_filename);
# $txt=extract_text_from_array(@searcharray);
# ($seq_query,$seq_subject)=get_seqs_from_queryline($query,$cons,$subject,$mask_nonconserved);

sub is_number{
    # returns true (1) if the input is a number
    # $t=is_number($record);

    my $record=$_[0];
    my $t=1;

    $record=~s/\.//;
    $record=~s/^\-//;
    if ($record =~ /\D/) {
	$t=0;
    }

    return $t;
}

sub extract_text_from_array{
    # usage:
    # ($txt,$n_matches)=extract_text_from_array($search_pattern,@searcharray);
    # where the search pattern is stored in $searcharray($#searchpattern);
    # if there is an equal sign, split it and return the number only

    # 09-14-2003: return n_matches as well
    # 07-16-2003: use search pattern directly as input (instead of inserting it into the array)

    my @searcharray;
    my $search_pattern;
    my @searchstuff;
    my @splitarray;
    my $txt="-1";
    my $temp;
    my $n_matches=0;

    ($search_pattern,@searcharray)=@_;

    @searchstuff=grep /$search_pattern/,@searcharray;
    if (@searchstuff) {
	$n_matches=$#searchstuff+1;
	if ($n_matches > 1) {
	    print "n_matches=$n_matches in extract_text_from_array\tparser_methods_v1.pl\n";
	}
	$temp=$searchstuff[0];
	if ($temp =~ /\=/) {
	    @splitarray=split /\=/,$searchstuff[0];
	    $txt=$splitarray[1];
	    if (length($txt)<1) {
		$txt="-1";
	    }
	} else {
	    $txt=$temp;
	}
    }
    chomp $txt;
    return ($txt,$n_matches);
}

sub extract_blast_results_from_array{
    # extract_blast_results_from_array($mask_nonconserved,@blast_output)
    # extract blast results from an array output of the form @blast_output=`blast blah blah blah`
    # usage: ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos)=extract_blast_results_from_array($mask_nonconserved,@blast_output);
    
    # 12-01-2003: added mask_nonconserved

    my $delimiter="Score";
    my $end_line;
    my $entry;
    my $expect="";
    my $i=0;
    my $identities_nt="";
    my $identities_perc="";
    my $n_lines;
    my $partial_expect;
    my $partial_identities_nt;
    my $partial_identities_perc;
    my $partial_qpos;
    my $partial_score;
    my $partial_seq_query;
    my $partial_seq_subject;
    my $partial_spos;
    my $query_pos="";
    my $score="";
    my $seq_query="";
    my $seq_subject="";
    my $start_line;
    my $subject_pos="";
    my @blast_output;
    my $mask_nonconserved;
    my @searchstuff;

    ($mask_nonconserved,@blast_output)=@_;

    $start_line=0;           # start from the very beginning, a good place to start
    $n_lines=$#blast_output; # number of lines in blast output
    $blast_output[$n_lines]=$delimiter;

    while ($start_line<$n_lines) {
	$blast_output[$n_lines+1]=$start_line;
	($entry,$end_line)=get_entry_beginning_with(@blast_output);

	$start_line=$end_line;
	($partial_score,$partial_expect,$partial_identities_nt,$partial_identities_perc,$partial_seq_query,$partial_seq_subject,$partial_qpos,$partial_spos)=extract_blast_results_from_entry($entry,0,$mask_nonconserved);
	#print "partial_score=$partial_score\n";
	#print "partial_expect=$partial_expect\n";
	#print "partial_identities_nt=$partial_identities_nt\n";
	#print "partial_identities_perc=$partial_identities_perc\n";
	#print "partial_seq=$partial_seq\n";
	if ($partial_identities_nt>0) {
	    $score.="$partial_score\n";
	    $expect.="$partial_expect\n";
	    $identities_nt.="$partial_identities_nt\n";
	    $identities_perc.="$partial_identities_perc\n";
	    $seq_query.="$partial_seq_query\n";
	    $seq_subject.="$partial_seq_subject\n";
	    $query_pos.="$partial_qpos\n";
	    $subject_pos.="$partial_spos\n";
	} else {
	    if ($i eq 0) {
		$identities_nt=0;
	    }
	}
	$i++;
    }
    #print "i=$i\n";
    return ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos);
}

sub extract_blastall_results_from_array{
    # extract_blastall_results_from_array($mask_nonconserved,@blast_output)
    # extract blastall results from an array output of the form @blast_output=`blast blah blah blah`
    # usage: ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos,$info)=extract_blastall_results_from_array($mask_nonconserved,@blast_output);
    # in constrat to extract_blast_results_from_array here we also extract the information line
    
    #my $delimiter="Score";
    my $delimiter1=">";
    my $delimiter2="Score";
    my $end_line;
    my $entry;
    my $expect="";
    my $i=0;
    my $identities_nt="";
    my $identities_perc="";
    my $n_lines;
    my $partial_expect;
    my $partial_identities_nt;
    my $partial_identities_perc;
    my $partial_qpos;
    my $partial_score;
    my $partial_seq_query;
    my $partial_seq_subject;
    my $partial_spos;
    my $query_pos="";
    my $score="";
    my $seq_query="";
    my $seq_subject="";
    my $start_line;
    my $subject_pos="";
    my @blast_output;
    my $mask_nonconserved;
    my @searchstuff;
    my @entry_array;
    my $entry_start;
    my $entry_n_lines;
    my $info_line;
    my $info="";

    ($mask_nonconserved,@blast_output)=@_;

    $start_line=0;           # start from the very beginning, a good place to start
    $n_lines=$#blast_output; # number of lines in blast output
    $blast_output[$n_lines]=$delimiter1;

    while ($start_line<$n_lines) {
	$blast_output[$n_lines+1]=$start_line;
	($entry,$end_line)=get_entry_beginning_with(@blast_output);
	$start_line=$end_line;
	
	@entry_array=split /\n/,$entry;
	$entry_n_lines=$#entry_array+1;
	for ($j=0;$j<$entry_n_lines;$j++) {
	    $entry_array[$j]="$entry_array[$j]\n";
	}
	$entry_start=0;
	$entry_array[$entry_n_lines]=$delimiter2;

	$info_line=$entry_array[0];
	$info_line=~s/^\>//;

	while ($entry_start<$entry_n_lines) {
	    $entry_array[$entry_n_lines+1]=$entry_start;
	    ($entry,$end_line)=get_entry_beginning_with(@entry_array);
	    $entry_start=$end_line;
	
	    ($partial_score,$partial_expect,$partial_identities_nt,$partial_identities_perc,$partial_seq_query,$partial_seq_subject,$partial_qpos,$partial_spos)=extract_blast_results_from_entry($entry,0,$mask_nonconserved);
	    #print "partial_score=$partial_score\n";
	    #print "partial_expect=$partial_expect\n";
	    #print "partial_identities_nt=$partial_identities_nt\n";
	    #print "partial_identities_perc=$partial_identities_perc\n";
	    #print "partial_seq=$partial_seq\n";
	    if ($partial_identities_nt>0) {
		$score.="$partial_score\n";
		$expect.="$partial_expect\n";
		$identities_nt.="$partial_identities_nt\n";
		$identities_perc.="$partial_identities_perc\n";
		$seq_query.="$partial_seq_query\n";
		$seq_subject.="$partial_seq_subject\n";
		$query_pos.="$partial_qpos\n";
		$subject_pos.="$partial_spos\n";
		$info.="$info_line";
	    } else {
		if ($i eq 0) {
		    $identities_nt=0;
		}
	    }
	    $i++;
	}
    }

    return ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos,$info);
}

sub get_entry_beginning_with{
    # get specific segment (e.g. a record) from a multi-line array input
    # this program derives from read_entry_beginning_with
    # the difference is that the latter uses a file as an input whereas this one uses a multi-line array input
    # usage: ($record,$end_line)=get_entry_beginning_with($delimiter,$start_line,@entry)

    #my $delimiter=$_[0];    # a piece of text which uniquely marks the end of a record (e.g. BLASTN   BLASTX   //   \\  etc...)
    #my $start_line=$_[1];   # initial processing line
    my $delimiter;
    my $start_line;
    my @entry=@_;           # input in a multi-line array format
    my $n_lines;
    my $end_line;
    my $input=0;
    my $record="";
    my $FIRST_TIME="TRUE";
    my $END_OF_RECORD="FALSE";
    my $current_file_position;
    my $curr_line;
   
    $n_lines=$#entry+1;

    $delimiter=$entry[$n_lines-2];
    $start_line=$entry[$n_lines-1];
    
    $n_lines=$n_lines-2;

    $curr_line=$start_line;

    while ( $END_OF_RECORD eq "FALSE"){
	$input=$entry[$curr_line];
	if ($curr_line>=$n_lines) {
	    #if (!($input)){ # end of file
	    $end_line=$n_lines;
	    return($record,$end_line);
	}
	if ($input =~ /$delimiter/) {
	    if ($FIRST_TIME eq "TRUE"){
		$FIRST_TIME = "FALSE";
	    } else {
		$END_OF_RECORD="TRUE";
		$end_line=$curr_line-1;
		return($record,$end_line); #send back what we got
	    }
	}
	if ($FIRST_TIME eq "FALSE") {
	    $record.=$input;
	}
	$curr_line++;
    } #endwhile        
    return($record,$end_line);
}#endsub

sub extract_blast_results_from_file{
    # extract_blast_results_from_file($blast_filename)
    # all blast results are reported
    # usage: ($score,$expect,$identities_nt,$identities_perc,$seq)=extract_blast_results_from_file($blast_filename);

    my $blast_file;
    my $blast_filename=$_[0];
    my $delimiter="Score";
    my $entry;
    my $expect="";
    my $i=0;
    my $identities_nt="";
    my $identities_perc="";
    my $partial_expect;
    my $partial_identities_nt;
    my $partial_identities_perc;
    my $partial_score;
    my $partial_seq;
    my $score="";
    my $seq="";
    my @searchstuff;

    if (-e $blast_filename) {
	$blast_file=open_read_file($blast_filename);
	while ($i<1000000) {
	    $entry=read_entry_beginning_with($blast_file,$delimiter);
	    if (! $entry) {
		#print "end of file\n";
		last;
	    }
	    #print "entry=$entry\n";
	    ($partial_score,$partial_expect,$partial_identities_nt,$partial_identities_perc,$partial_seq)=extract_blast_results_from_entry($entry,0,0);
	    #print "partial_score=$partial_score\n";
	    #print "partial_expect=$partial_expect\n";
	    #print "partial_identities_nt=$partial_identities_nt\n";
	    #print "partial_identities_perc=$partial_identities_perc\n";
	    #print "partial_seq=$partial_seq\n";
	    if ($partial_identities_nt>0) {
		$score.="$partial_score\n";
		$expect.="$partial_expect\n";
		$identities_nt.="$partial_identities_nt\n";
		$identities_perc.="$partial_identities_perc\n";
		$seq.="$partial_seq\n";
	    } else {
		if ($i eq 0) {
		    $identities_nt=0;
		}
	    }
	    $i++;
	}
    }
    #print "i=$i\n";
    return ($score,$expect,$identities_nt,$identities_perc,$seq);
}

sub extract_blast_results_from_entry{
    # extract blast results (score, expectation, identities, strand, positions) from entry of the form:
    #
    # Score =  107 bits (54), Expect = 1e-28
    # Identities = 54/54 (100%)
    # Strand = Plus / Plus
    #
    #                                                            
    # Query: 4  gctactactttgaagctgttgcgcagctgccgcaggagacgcgcaaccagctgg 57
    #           ||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # Sbjct: 32 gctactactttgaagctgttgcgcagctgccgcaggagacgcgcaaccagctgg 85
    #
    # usage:
    # ($score,$expect,$identities_nt,$identities_perc,$seq_query,$seq_subject,$query_pos,$subject_pos)=extract_blast_results_from_entry($entry,$v,$mask_nonconserved);

    # 12-01-2003: add mask_nonconserved, if this is 1, then return "n" in the non-conserved places for the sequences
    # 02-06-2003: add query and subject positions as output
    # 06-19-2003: return seq_query and seq_subject

    my $entry=$_[0];
    my @searchstuff;
    my @subjectsearch;
    my $score="na";
    my $expect="na";
    my $identities_nt=0;
    my $identities_perc=0;
    my $seq1="na";
    my $seq2="na";
    my $query_pos="na";
    my $subject_pos="na";
    my $partial_seq1;
    my $partial_seq2;
    my $queryline;
    my $n_querylines;
    my $i;
    my @entry_array;
    my $v=$_[1]; # verbose or not
    my $temp_pos;
    my $q;
    my $c;
    my $s;
    my $record;
    my $mask_nonconserved=$_[2];

    # read info 
    (@searchstuff)=$entry=~ /Score\s=\s+(\w+\.?\w?) bits/;
    if (!@searchstuff) {
	#print "error! i could not find score\nentry=$entry\n";
	if ($v eq 1) {
	    print "warning! i could not find score\n";
	} 
    } else {
	$score=$searchstuff[0];
	(@searchstuff)=$entry=~/Expect = (\w+\.?\-?\w+)/g;
	if (!@searchstuff) {
	    if ($v eq 1) {
		print "error! i could not find expected value\nentry=$entry\n";
	    }
	} else {
	    $expect=$searchstuff[0];
	    if ($expect =~ m/^e/ ) {
		$expect="1$expect";
	    }
	    (@searchstuff)=$entry=~ /Identities = (\w+)\//;
	    if (!@searchstuff) {
		print "warning! i could not find the Identities line\n";
	    } else {
		$identities_nt=$searchstuff[0];
	    }		
	    (@searchstuff)=$entry=~ /\((\w+)\%\)/;
	    if (!@searchstuff) {
		print "warning! i could not find the Identities percentage information\n";
	    } else {
		$identities_perc=$searchstuff[0];
	    }
	    # read sequence
	    $seq1="";
	    $seq2="";
	    $query_pos="";
	    $subject_pos="";
	    @entry_array=split /\n/,$entry;
	    (@searchstuff)=grep /Query/,@entry_array;
	    (@subjectsearch)=grep /Sbjct/,@entry_array;
	    #(@searchstuff)=$entry=~ /Query: \d+ (\w+) /;
	    if (!@searchstuff) {
		print "error! i could not find query line\n";
	    } else {
		$n_querylines=$#searchstuff+1;

		# get initial query and subject positions
		$queryline=$searchstuff[0];
		$temp_pos=get_pos_from_queryline($queryline,1);
		$query_pos="${temp_pos}";
		$queryline=$subjectsearch[0];
		$temp_pos=get_pos_from_queryline($queryline,2);
		$subject_pos="${temp_pos}";

		# get final query and subject positions
		$queryline=$searchstuff[$n_querylines-1];
		$temp_pos=get_pos_from_queryline($queryline,3);
		$query_pos="${query_pos}\t${temp_pos}";
		$queryline=$subjectsearch[$n_querylines-1];
		$temp_pos=get_pos_from_queryline($queryline,3);
		$subject_pos="${subject_pos}\t${temp_pos}";
	    }

	    # get sequences
	    for ($i=0;$i<$#entry_array;$i++) {
		$record=$entry_array[$i];
		if ($record =~ /^Query/) {
		    $q=$entry_array[$i];            # query sequence
		    $c=$entry_array[$i+1];          # conservation line
		    $s=$entry_array[$i+2];          # subject sequence
		    
		    ($partial_seq1,$partial_seq2)=get_seqs_from_queryline($q,$c,$s,$mask_nonconserved);
		    $seq1.=$partial_seq1;
		    $seq2.=$partial_seq2;
		}
	    }
	}
    }

    return ($score,$expect,$identities_nt,$identities_perc,$seq1,$seq2,$query_pos,$subject_pos);
}

sub get_seqs_from_queryline{
    # get both subject and query sequences
    # if mask_nonconserved==1, then return only conserved segments, "n" masking the rest
    # use only a single line of query entries
    # usage:
    # ($seq_query,$seq_subject)=get_seqs_from_queryline($query,$cons,$subject,$mask_nonconserved);
    #
    # 12-01-2003: add mask_nonconserved, if this is 1 then replace non-conserved nucleotides by "n"

    my $query=$_[0];
    my $cons_orig=$_[1];
    my $cons=$cons_orig;
    my $subject=$_[2];
    my $mask_nonconserved=$_[3];
    my @query_array;
    my @subject_array;
    my @cons_array;
    my $i;
    my $n;
    my $seq1=get_sequence_from_queryline($query);
    my $seq2=get_sequence_from_queryline($subject);
    my $seq_query="";
    my $seq_subject="";
    my $n1=length($query);
    my $temp=$query;
    my $n2;
    my $n3;

    $temp=~s/Query\: //;
    $temp=~s/\d+//;
    $temp=~s/\s+//;
    $n2=length($temp);
    $n3=$n1-$n2;

    for ($i=0;$i<$n3;$i++) {
	$cons=~s/\s//;
    }

    (@query_array)=string2array($seq1);
    (@subject_array)=string2array($seq2);
    (@cons_array)=string2array($cons);

    $n=$#query_array;
    for ($i=0;$i<$n;$i++) {
	if ($masked_nonconserved eq 1) {
	    if ($cons_array[$i] eq "|") {
		$seq_query.="$query_array[$i]";
		$seq_subject.="$subject_array[$i]";
		if ($query_array[$i] eq "-") {
		    print "WARNING!\n";
		    print "$i\t$query_array[$i]\n";
		    print "seq1\n$seq1\n";
		    print "seq2\n$seq2\n";
		    print "cons\n$cons\n";
		    print "query\n$query\n";
		    print "subject\n$subject\n";
		    print "cons_orig\n$cons_orig\n";
		    print "n1=$n1\tn2=$n2\tn3=$n3\n";
		    exit;
		}
	    } else {
		if ($query_array[$i] ne "-") { 
		    $seq_query.="n";
		}
		if ($subject_array[$i] ne "-") {
		    $seq_subject.="n";
		}
	    }
	} else {
	    if ($query_array[$i] ne "-") {
		$seq_query.="$query_array[$i]";
	    }
	    if ($subject_array[$i] ne "-") {
		$seq_subject.="$subject_array[$i]";
	    }
	}
    }

    return ($seq_query,$seq_subject);
}

sub get_pos_from_queryline{
    # get position from query or subject line
    # usage: $pos=get_pos_from_queryline($queryline,$querysubj);
    # where $querysubj=1 for query, 2 for subject and 3 for ending position
    # e.g.
    # Query: 4  gctactactttgaagctgttgcgcagctgccgcaggagacgcgcaaccagctgg 57

    my $queryline=$_[0];
    my $querysubj=$_[1];
    my @searchstuff;
    my $searchpat;
    my $pos=-1;

    if ($querysubj eq 1) {
	$searchpat='Query\:\s+(\d+)\s+';
    } else {
	if ($querysubj eq 2) {
	    $searchpat='Sbjct\:\s+(\d+)\s+';
	} else {
	    if ($querysubj eq 3) {
		$searchpat='[a|c|g|t|n]\s+(\d+)';
	    }
	}
    }

    (@searchstuff) = $queryline =~ /${searchpat}/;
    if (@searchstuff) {
	$pos=$searchstuff[0];
    }
    
    return $pos;
}


sub get_sequence_from_queryline{
    # get sequence from a blast query line in the form:
    # Query: 4  gctactactttgaagctgttgcgcagctgccgcaggagacgcgcaaccagctgg 57
    # 
    # usage:
    # $seq=get_sequence_from_queryline($queryline);

    my $queryline=$_[0];
    my $seq;
    my @splitarray;

    @splitarray=split /\s+/,$queryline;
    $seq=$splitarray[2];

    return ($seq);
}

sub split_infowm{
    # split a fasta formatted weight matrix entry as in 
    # > info about the motif
    # frequency matrix for a,c,g,t
    # into an info variable with the first line
    # and a wm variable with the next 4 lines
    # return also the motif length
    # usage:
    # ($info,$wm,$ml)=split_infowm($entry);

    my $entry=$_[0];
    my @splitarray=split /\n/,$entry;
    my $info=$splitarray[0];
    my $wm="";
    my $j;
    my $ml;

    for ($j=1;$j<=4;$j++) {
	$wm.="$splitarray[$j]\n";
    }	
    $entry=$splitarray[1];
    @splitarray=split /\t/,$entry;
    $ml=$#splitarray+1;

    return ($info,$wm,$ml);
}

sub split_infoseq{
    # split a fasta formated entry as in 
    # >info about the transcript
    # accccaca
    # into an info variable with the first line
    # and a seq variable with the sequence (in a single line)
    # usage: ($info,$seq)=split_infoseq($entry);

    my $entry=$_[0];
    my $info;
    my $seq;
    my @splitarray;
    my $n;
    my $i;

    $seq="";
    @splitarray=split /\n/,$entry;
    $n=$#splitarray+1;
    $info=$splitarray[0];
    for ($i=1;$i<$n;$i++) {
	$seq.=$splitarray[$i];
    }

    return ($info,$seq);
}

sub read_entry_beginning_with{
    # this was written by keith ching
    #       function: to read in a "record" from a flat file like genbank,
    #       unigene, pfam, swissprot, etc and return that record as a scalar

    my $filehandle=$_[0];   # of input file handle to search for entry
    my $delimiter=$_[1];    # a piece of text which uniquely marks the end of a file
                                # ie BLASTN   BLASTX   //   \\  etc...
    my $input=0;
    my $record;
    my $FIRST_TIME="TRUE";
    my $END_OF_ENTRY="FALSE";
    my $current_file_position;
    while ( $END_OF_ENTRY eq "FALSE"){
	$input=<$filehandle>;
	if (!($input)){ # end of file
	    return($record);
	}
	if ($input=~/$delimiter/){
	    if ($FIRST_TIME eq "TRUE"){
		$FIRST_TIME = "FALSE";
	    }else{
		$END_OF_ENTRY="TRUE";
		seek $filehandle, $current_file_position, 0; # reset the file position to the last line
		return($record); #send back what we got
	    }
	}
	$record.=$input;
	$current_file_position=tell $filehandle; # save the position of the last line
    } #endwhile        
    return($record);
}#endsub

sub gk_read_entry_beginning_with{
    # i modified the previous version to skip any text at the beginning of the file before the first entry
    #       function: to read in a "record" from a flat file like genbank,
    #       unigene, pfam, swissprot, etc and return that record as a scalar

    my $filehandle=$_[0];   # of input file handle to search for entry
    my $delimiter=$_[1];    # a piece of text which uniquely marks the end of a file
                                # ie BLASTN   BLASTX   //   \\  etc...
    my $input=0;
    my $record;
    my $FIRST_TIME="TRUE";
    my $END_OF_ENTRY="FALSE";
    my $current_file_position;
    while ( $END_OF_ENTRY eq "FALSE"){
	$input=<$filehandle>;
	if (!($input)){ # end of file
	    return($record);
	}
	if ($input=~/$delimiter/){
	    if ($FIRST_TIME eq "TRUE"){
		$FIRST_TIME = "FALSE";
		$record="";
	    }else{
		$END_OF_ENTRY="TRUE";
		seek $filehandle, $current_file_position, 0; # reset the file position to the last line
		return($record); #send back what we got
	    }
	}
	$record.=$input;
	$current_file_position=tell $filehandle; # save the position of the last line
    } #endwhile        
    return($record);
}#endsub

sub search_pattern{
    my $entry=$_[0];                #text to search
    my $searchpattern=$_[1];
    my @query;
    
    @query=$entry=~/$searchpattern/g;
    return(@query);
}#endsub

sub singleline2multipleline{
    # add line breaks to convert a single line to multiple lines
    # usage: $txt2=singleline2multipleline($txt1,$n);

    my $inputtxt=$_[0];
    my $n=$_[1];
    my $outputtxt="";
    my $i;

    for ($i=0;$i<=length($inputtxt)-$n;$i=$i+$n) {
	$outputtxt.=substr($inputtxt,$i,$n);
	$outputtxt.="\n";
    }
    $outputtxt.=substr($inputtxt,$i,length($inputtxt)%$n);

    return ($outputtxt);
}

1;   #it took me about one hour to realise i had to put this guy here!

