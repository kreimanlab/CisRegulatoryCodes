#!/usr/bin/perl
#fileio_methods_v1.pl
#methods (subs) used to deal with files
#
# 11-05-2003: moved is_comment from regul_search_methods.pl to here
# 11-05-2003: corrected a bug in read_whole_file_v4 so as to accept comment lines in the middle or end of the file as well

# list of methods
# ($system_status,$warning_msg)=my_rename($filename_orig,$filename_new,$verbose,$exit_on_error);
# $FILEHANDLE=open_read_file($filename);
# @filecontents=read_whole_file($filename);
# @filecontents=read_whole_file_v2($filename,$header);
# @filecontents=read_whole_file_v3($filename,$header,$addline);
# @filecontents=read_whole_file_v4($filename,$header,$addline);
# $sec=convert_time_to_sec($day1,$day2,$hour1,$hour2,$min2,$min1,$sec1,$sec2); 
# $filename_noext=get_filename_noext($filename_withext);
# @columncontents=read_column($filename,$column,$header,$addline);
# @filecontents=read_partial_file($filename,$symbol,$header,$addline);
# $linenr=get_line_number($textline);
# $coltext=get_column($textline,$col);
# ($n,$n_comments)=line_count($filename);
# sort_file_num($input_filename,$column);
# $status=copyfile2handle($file2copy,$filehandle,$copycomments);
# @columncontents=read_column_v2($filename,$column,$header,$addline);
# @columncontents=read_column_v3($filename,$column,$header,$addline,$separator);
# $header=is_header($txt);
# $t=get_time_string();
# @fasta=read_whole_fasta($filename);
# $n=fasta_line_count($filename);
# @matrix=read_matrix_file($filename,$n_cols);
# ($entry,$info,$errorcode)=extract_seq_from_fasta($filename,$id,$returnfasta);
# $truefalse=is_comment($record);

sub gk{
    # gk
    print "\ngabriel kreiman\n";
    print "kreiman@mit.edu\n";
    print "all rights reserved (c)\n";
}

sub separate_dir_name
{
    # given a file name, return the directory and the file name without the directory
    # usage:
    # ($folder,$just_filename)=separate_dir_name($filename);

    my $filename=$_[0];
    my $ri;
    my $just_filename;
    my $folder;
    my $l;

    $ri=rindex($filename,"/");
    if ($ri > 0) {
	$l=length($filename);
	$just_filename=substr($filename,$ri+1,$l);
	$folder=substr($filename,0,$ri+1);
    } else {
	$just_filename=$filename;
	$folder="";
    }
    
    return ($folder,$just_filename);
}

sub get_file_modifiedt{
    # given a file, return its modified date
    
    my $filename=$_[0];
    my $month;
    my $day;
    my $year;
    my @splitarray;
    my @loutput;
    my $tstring;

    if (-e $filename) {
	$code_exec="ls -l $filename";
	@loutput=`$code_exec`;
	@splitarray=split /\s+/,$loutput[0];
	$month=$splitarray[5];
	$month=month2number($month);
	$day=$splitarray[6];
	$year=$splitarray[7];
	if ($year =~ /\:/ ) {
	    $year=2004;
	}
    } else {
	#print "$filename does not exist\n";
	$month=-1;
	$day=-1;
	$year=-1;
    }
    $tstring="${month}${day}${year}";

    return ($year,$month,$day,$tstring);
}

sub month2number{
    # convert a string month to a number
    my $month_string=$_[0];
    my $month_number=-1;

    if ($month_string eq "Jan") {
	$month_number=1;
    }
    if ($month_string eq "Feb") {
	$month_number=2;
    }
    if ($month_string eq "Mar") {
	$month_number=3;
    }
    if ($month_string eq "Apr") {
	$month_number=4;
    }
    if ($month_string eq "May") {
	$month_number=5;
    }
    if ($month_string eq "Jun") {
	$month_number=6;
    }
    if ($month_string eq "Jul") {
	$month_number=7;
    }
    if ($month_string eq "Aug") {
	$month_number=8;
    }
    if ($month_string eq "Sep") {
	$month_number=9;
    }
    if ($month_string eq "Oct") {
	$month_number=10;
    }
    if ($month_string eq "Nov") {
	$month_number=11;
    }
    if ($month_string eq "Dec") {
	$month_number=12;
    }
    return $month_number;
}

sub is_comment{
    # $truefalse=is_comment($record);
    # returns "true" if record starts with ">" or %

    my $record=$_[0];
    my $truefalse="";

    if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
	$truefalse="true";
    } else {
	$truefalse="false";
    }
    return $truefalse;
}

sub my_rename
{
    # run rename command
    # usage:
    # ($system_status,$warning_msg)=my_rename($filename_orig,$filename_new,$verbose,$exit_on_error);

    my $filename_orig=$_[0];
    my $filename_new=$_[1];
    my $verbose=$_[2];
    my $exit_on_error=$_[3];
    my $system_status=0;
    my $warning_msg="";
    my $n_renamed;
    
    if (-e $filename_new) {
	$warning_msg.="WARNING! $filename_new already exists. overwriting...\n";
    }
    $n_renamed=rename($filename_orig,$filename_new);
    if ($n_renamed ne 1) {
	if ($exit_on_error eq 1) {
	    print "error!\nfilename_orig=$filename_orig\nfilename_new=$filename_new\nn_renamed=$n_renamed\n";
	    exit;
	} else {
	    $warning_msg.="error!\tfilename_orig=$filename_orig\tfilename_new=$filename_new\tn_renamed=$n_renamed\n";
	    $system_status=-1;
	    if ($verbose) {
		print "$warning_msg\n";
	    }
	}
    } else {
	$warning_msg.="\tmoved $filename_orig to $filename_new OK\n";
	if ($verbose) {
	    print "$warning_msg\n";
	}
    }

    return ($system_status,$warning_msg);
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

sub read_matrix_file{
    # read a file with a matrix and store each column separately
    # usage:
    # @matrix=read_matrix_file($filename,$n_cols);

    my $filename=$_[0];
    my $n_cols=$_[1];
    my $FILEHANDLE;
    my $record;
    my @matrix;
    my $i=0;
    my $j;
    my $read_header;
    my @splitarray;

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    $read_header=1;
    while ($read_header eq 1) {
	$record=<$FILEHANDLE>;
	if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
	} else {
	    $read_header=0;
	    chomp $record;
	    $i++;
	    @splitarray=split /\t/,$record;
	    for ($j=0;$j<$n_cols;$j++) {
		$matrix[$i][$j+1]=$splitarray[$j];
	    }
	}
    }

    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	@splitarray=split /\t/,$record;
	for ($j=0;$j<$n_cols;$j++) {
	    $matrix[$i][$j+1]=$splitarray[$j];
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @matrix;
}

sub get_filenumber{
    # usage:
    # $n=get_filenumber($prefix,$suffix);

    my $prefix=$_[0];
    my $suffix=$_[1];
    my $i=1;
    my $gotit=1;
    my $filename;

    while ($gotit eq 1) {
	$filename="${prefix}${i}${suffix}";
	if (-e $filename) {
	    $i++;
	} else {
	    $gotit=0;
	}
    }

    return $i;
}

sub get_time_string{
    # returns a pretty time string for log purposes
    # usage: $t=get_time_string();
    my $Second;
    my $Minute;
    my $Hour;
    my $Year;
    my $WeekDay;
    my $DayOfYear;
    my $IsDST;
    my $t;
    
    ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
    $Year+=1900;
    $Month++;

    $t="$Month \/ $Day \/ $Year $Hour\:$Minute\:$Second";

    return $t;
}

sub read_whole_fasta{
    # read all sequences from a fasta file
    # usage: @fasta=read_whole_fasta($filename);
    # returns all the fasta file such that $fasta[$i][1]=info and $fasta[$i][2]=sequence (in fasta format) for the i-th entry
    
    my $input_filename=$_[0];
    my $input_file=open_read_file($input_filename);
    my $entry;
    my $i=0;
    my $info;
    my $sequence;
    my @fasta;

    $fasta[0][1]=$input_filename;
    while ($i<1000000) {
	$entry=read_entry_beginning_with($input_file,">");
	if (!$entry) {
	    last;
	}
	($info,$sequence)=split_infoseq($entry);
	$i++;
	$fasta[$i][1]=$info;
	$sequence=singleline2multipleline($sequence,80);
	$fasta[$i][2]=$sequence;
    }
    close ($input_file);

    return @fasta;
}

sub is_header{
    # usage:
    # $header=is_header($txt);
    # returns 1 if $txt starts with ">" or "%"

    my $record=$_[0];
    my $header=0;

    if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
	$header=1;
    }
    return $header;
}

sub copyfile2handle{
    # copy contents of a file to another one given its file handle
    # usage:
    # $status=copyfile2handle($file2copy,$filehandle,$copycomments);
    # copy comment lines also if $copycomments is 1
    #
    # 09-29-2003: added close(input_file);

    my $file2copy=$_[0];
    my $filehandle=$_[1];
    my $copycomments=$_[2];
    my $status=0;
    my $record;

    if (-e $file2copy) {
	open (input_file,$file2copy) || die "could not open $file2copy";
	if ($copycomments eq 1) {
	    while ($record=<input_file>) {
		print $filehandle "$record";
	    }
	} else {
	    while ($record=<input_file>) {
		if ( ($record =~ /^\>) | ($record =~ /^\%/) ) {
		} else {
		    print $filehandle "$record";
		}
	    }
	}
	close (input_file);
    } else {
	$status=-1;
    }
    
    return ($status);
}

sub sort_file_num{
    # sort a file by a numeric column
    # usage:
    # sort_file_num($input_filename,$column);
    #
    # 11-05-2003: changed <=n_rows to <n_rows to avoid printing an empty line at the end

    my $input_filename=$_[0];
    my $input_filename_noext;
    my $output_filename;
    my $col=$_[1];
    my $i;
    my $j;
    my $record;
    my $n_cols;
    my $n_rows;
    my @rank;
    my $curr_index;
    my @data;
    my @data_col;
    my @splitarray;
    my @data_sorted;
 
    $input_filename_noext=get_filename_noext($input_filename);
    $output_filename=">${input_filename_noext}.sortedbycol${col}.txt";

    $col--;
    
    open (input_file,$input_filename) || die "could not open ${input_filename}";
    open (sortfilenum_output,$output_filename) || die "could not open ${output_filename}";

    while ($record=<input_file>) {
	chomp $record;
	if ( ($record =~ /^\>/)  | ($record =~ /^\%/) ) {
	    print sortfilenum_output "$record\n";
	} else {
	    @splitarray=split /t/,$record;
	    $n_cols=$#splitarray;
	    for ($j=0;$j<=$#splitarray;$j++) {
		$data[$i][$j]=$splitarray[$j];
	    }
	    $data_col[$i]=$splitarray[$col];
	    $i++;
	}
    }
    close (input_file);
    $n_rows=$i;

    @rank[sort { $data_col[$a] <=> $data_col[$b] } 0 .. $#data_col] = 0 .. $#data_col;

    for ($i=0;$i<$n_rows;$i++) {
	$curr_index=$rank[$i];
	for ($j=0;$j<=$n_cols;$j++) {
	    $data_sorted[$curr_index][$j]=$data[$i][$j];
	}
    }

    #for ($i=0;$i<=$n_rows;$i++) {
    for ($i=0;$i<$n_rows;$i++) {
	for ($j=0;$j<=${n_cols};$j++) {
	    print sortfilenum_output "$data_sorted[$i][$j]\t";
	}
	print sortfilenum_output "\n";
    }
    close (sortfilenum_output);
}

sub fasta_line_count{
    # count number of entries in a fasta-formatted file
    # usage:
    # $n=fasta_line_count($filename);
    #
    
    my $filename=$_[0];
    my $n=0;
    my $entry;
    my $file;
    
    $file=open_read_file($filename);
    while ($n<1000000) {
	$entry=read_entry_beginning_with($file,">");
	if (!$entry) {
	    last;
	}
	$n++;
    }
    close ($file);

    return $n;
}

sub line_count{
    # count the number of lines in a file
    # ($n,$n_comments)=line_count($filename);
    my $filename=$_[0];
    my $n=0;
    my $record;
    my $n_comments=0;
    
    if (-e $filename) {
	open (input_file,$filename) || die "could not open $filename for reading";
	while ($record=<input_file>) {
	    if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
		$n_comments++;
	    }
	    $n++;
	}
    }
    return ($n,$n_comments);
}

sub conservative_rename{
    # if the new file name exists, convert it to old before renaming to avoid overwriting
    # usage:
    # conservative_rename($f_original,$f_target);

    my $f1=$_[0];
    my $f2=$_[1];
    my $f_temp;
    my $n_renamed;

    if (-e $f2) {
	$f_temp="$f2.old";
	print "moving $f2 to $f_temp...";
	$n_renamed=rename($f2,$f_temp);
	if ($n_renamed eq 1) {
	    print "ok\n";
	} else {
	    print "error\n";
	}
    }

    $n_renamed=rename($f1,$f2);
    if ($n_renamed eq 1) {
	print "moving $f1 to $f2... ok\n";
    } else {
	print "moving $f1 to $f2... error\n";
    }
}

sub get_line_number{
    # get line number (first entry) from a tab separated text
    # usage:
    # $linenr=get_line_number($textline);

    my $textline=$_[0];
    my @splitarray;
    my $linenr;

    @splitarray=split /\t/,$textline;
    $linenr=$splitarray[0];

    return $linenr;
}

sub get_column{
    # get column from a tab separated text line
    # usage:
    # $coltext=get_column($textline,$col);

    my $textline=$_[0];
    my $col=$_[1];
    my @splitarray;
    my $coltext;

    @splitarray=split /\t/,$textline;
    $coltext=$splitarray[$col];

    return $coltext;
}

sub get_filename_noext{
    # get_filename_noext
    # usage: $filename_noext=get_filename_noext($filename_withext);
    # return the name of the file without the extension (defined as the last 4 characters)

    # 11-08-2003 use rindex to locate the dot
    # 05-07-2003 incorporate fasta sequences
			   			   
    my $filename_wext=$_[0];
    my $filename_noext;
    my $l;
    my $ri;
    
    if (-e $filename_wext) {
	$ri=rindex($filename_wext,".");
	if ($r1<0) {
	    $filename_noext=$filename_wext;
	} else {	
	    $filename_noext=substr($filename_wext,0,$ri);
	} 
    } else {
	$filename_noext="none";
    }

    return ($filename_noext);
}

sub open_read_file{
    my $filename=$_[0];
    my $FILEHANDLE;
    my $temp;
    
    #open ($FILEHANDLE,$filename) || die "could not open file: $filename : $!";
    open ($FILEHANDLE,$filename) || die "could not open file: $filename; exiting...";

    return $FILEHANDLE;
}#endsub

sub read_whole_file{
    my $filename=$_[0];
    my $FILEHANDLE;
    my $record;
    my @filecontents;
    my $i=0;

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    $filecontents[0]=$filename;
    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	$filecontents[$i]=$record;
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @filecontents;
}

sub read_whole_file_v2{
    # new in v2, process header line separately if it exists
    # usage: @filecontents=read_whole_file_v2($filename,$header);
    # where $header = 1 if there is a header line
    my $filename=$_[0];
    my $header=$_[1];
    my $FILEHANDLE;
    my $record;
    my @filecontents;
    my $i=0;

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    if ($header eq 0) {
	$filecontents[0]=$filename;
    } else {
	$filecontents[0]=<$FILEHANDLE>;
    }
    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	$filecontents[$i]=$record;
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @filecontents;
}

sub read_whole_file_v3{
    # new in v3, optionally add line number
    # new in v2, process header line separately if it exists
    # usage: @filecontents=read_whole_file_v3($filename,$header,$addline);
    # where $header = 1 if there is a header line
    my $filename=$_[0];
    my $header=$_[1];
    my $addline=$_[2];
    my $FILEHANDLE;
    my $record;
    my @filecontents;
    my $i=0;

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    if ($header eq 0) {
	$filecontents[0]=$filename;
    } else {
	if ($addline eq 1) {
	    $i++;
	    $filecontents[0]="$i\t<$FILEHANDLE>";
	} else {
	    $filecontents[0]=<$FILEHANDLE>;
	}
    }
    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	if ($addline eq 1) {
	    $filecontents[$i]="$i\t$record";
	} else {
	    $filecontents[$i]=$record;
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @filecontents;
}

sub read_whole_file_v4{
    # new in v4, ignore all header lines, not just one; also allow matlab header lines
    # new in v3, optionally add line number
    # new in v2, process header line separately if it exists
    # usage: @filecontents=read_whole_file_v4($filename,$header,$addline);
    # where $header = 1 if there is a header line
    #
    # 11-05-2003 modify to allow comment lines after the meat of the file

    my $filename=$_[0];
    my $header=$_[1];
    my $addline=$_[2];
    my $FILEHANDLE;
    my $record;
    my @filecontents;
    my $i=0;
    my $read_header;
    my $temp_header="";

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    if ($header eq 0) {
	$filecontents[0]=$filename;
    } else {
	$read_header=1;
	while ($read_header eq 1) {
	    $record=<$FILEHANDLE>;
	    if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
		$temp_header.="$record";
	    } else {
		$read_header=0;
		#$filecontents[0]=$temp_header;
		chomp $record;
		$i++;
		if ($addline eq 1) {
		    $filecontents[$i]="$i\t$record";
		} else {
		    $filecontents[$i]="$record";
		}
	    }
	}
    }

    while ($record=<$FILEHANDLE>) {
	chomp $record;
	chomp $record;
	if ( ($record =~ /^\>/) | ($record =~ /^\%/) ) {
	    $temp_header.="$record\n";
	} else {
	    $i++;
	    if ($addline eq 1) {
		$filecontents[$i]="$i\t$record";
	    } else {
		$filecontents[$i]="$record";
	    }
	}
    }
    if ($header eq 1) {
	$filecontents[0]=$temp_header;
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @filecontents;
}

sub read_partial_file{
    # read selected entries (starting with a given symbol) from a file
    # usage: @filecontents=read_partial_file($filename,$symbol,$header,$addline);
    # where $header = 1 if there is a header line
    # $addline = 1 to add line number infomration
    my $filename=$_[0];
    my $symbol=$_[1];
    my $header=$_[2];
    my $addline=$_[3];
    my $FILEHANDLE;
    my $record;
    my @filecontents;
    my $i=0;

    $FILEHANDLE=open_read_file($filename);
    print "reading $filename ...";
    if ($header eq 0) {
	$filecontents[0]=$filename;
    } else {
	if ($addline eq 1) {
	    #$i++;
	    $filecontents[0]="$i\t<$FILEHANDLE>";
	} else {
	    $filecontents[0]=<$FILEHANDLE>;
	}
    }
    while ($record=<$FILEHANDLE>) {
	if ($record =~ /^${symbol}/) {
	    chomp $record;
	    $i++;
	    if ($addline eq 1) {
		$filecontents[$i]="$i\t$record";
	    } else {
		$filecontents[$i]=$record;
	    }
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @filecontents;
}

sub read_column{
    # read a column from a file
    # usage: @columncontents=read_column($filename,$column,$header,$addline);
    # where 
    # $header = 1 if there is a header line
    # $addline = 1 to add the line number
    # $column is the column number (we subtract one in here)
    # here we assume tab separated columns
    my $filename=$_[0];
    my $column=$_[1];
    my $header=$_[2];
    my $addline=$_[3];
    my $FILEHANDLE;
    my $record;
    my @columncontents;
    my $i=0;
    my @splitarray;

    $column--; # to account for the fact that it starts at 0
    $FILEHANDLE=open_read_file($filename);
    print "reading column $column from $filename ...";
    if ($header eq 0) {
	$columncontents[0]=$filename;
    } else {
	$record=<$FILEHANDLE>;
	if ($addline eq 1) {
	    #$i++;
	    $columncontents[0]="$i\t$record";
	} else {
	    $columncontents[0]=$record;
	}
    }
    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	@splitarray=split /\t/,$record;
	if ($addline eq 1) {
	    $columncontents[$i]="$i\t$splitarray[$column]";
	} else {
	    $columncontents[$i]=$splitarray[$column];
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @columncontents;
}

sub read_column_v2{
    # read a column from a file
    # usage: @columncontents=read_column_v2($filename,$column,$header,$addline);
    # where 
    # $header = 1 to skip any header lines
    # $addline = 1 to add the line number
    # $column is the column number (we subtract one in here)
    # here we assume tab separated columns

    my $filename=$_[0];
    my $column=$_[1];
    my $header=$_[2];
    my $addline=$_[3];
    my $FILEHANDLE;
    my $record;
    my @columncontents;
    my $i=0;
    my @splitarray;
    my $read_header;
    my $temp_header="";

    $column--; # to account for the fact that it starts at 0
    $FILEHANDLE=open_read_file($filename);
    print "reading column $column from $filename ...";
    if ($header eq 0) {
	$columncontents[0]=$filename;
    } else {
	$read_header=1;
	while ($read_header eq 1) {
	    $record=<$FILEHANDLE>;
	    chomp $record;
	    if (is_header($record)) {
		$temp_header.="$record\t";
	    } else {
		$read_header=0;
		@splitarray=split /\t/,$record;
		if ($addline eq 1) {
		    $columncontents[0]="$i\t$temp_header";
		    $i++;
		    $columncontents[1]="$i\t$splitarray[$column]";
		} else {
		    $columncontents[0]=$temp_header;
		    $i++;
		    $columncontents[1]=$splitarray[$column];
		}
		
	    }
	}
    }

    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	@splitarray=split /\t/,$record;
	if ($addline eq 1) {
	    $columncontents[$i]="$i\t$splitarray[$column]";
	} else {
	    $columncontents[$i]=$splitarray[$column];
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @columncontents;
}

sub read_column_v3{
    # read a column from a file
    # usage: @columncontents=read_column_v2($filename,$column,$header,$addline,$separator);
    # where 
    # $header = 1 to skip any header lines
    # $addline = 1 to add the line number
    # $column is the column number (we subtract one in here)
    # here we assume tab separated columns
    # new in v3: use an arbitrary separator

    my $filename=$_[0];
    my $column=$_[1];
    my $header=$_[2];
    my $addline=$_[3];
    my $separator=$_[4];
    my $FILEHANDLE;
    my $record;
    my @columncontents;
    my $i=0;
    my @splitarray;
    my $read_header;
    my $temp_header="";

    $column--; # to account for the fact that it starts at 0
    $FILEHANDLE=open_read_file($filename);
    print "reading column $column from $filename ...";
    if ($header eq 0) {
	$columncontents[0]=$filename;
    } else {
	$read_header=1;
	while ($read_header eq 1) {
	    $record=<$FILEHANDLE>;
	    chomp $record;
	    if (is_header($record)) {
		$temp_header.="$record\t";
	    } else {
		$read_header=0;
		@splitarray=split /${separator}/,$record;
		if ($addline eq 1) {
		    $columncontents[0]="$i\t$temp_header";
		    $i++;
		    $columncontents[1]="$i\t$splitarray[$column]";
		} else {
		    $columncontents[0]=$temp_header;
		    $i++;
		    $columncontents[1]=$splitarray[$column];
		}
		
	    }
	}
    }

    while ($record=<$FILEHANDLE>) {
	chomp $record;
	$i++;
	@splitarray=split /${separator}/,$record;
	if ($addline eq 1) {
	    $columncontents[$i]="$i\t$splitarray[$column]";
	} else {
	    $columncontents[$i]=$splitarray[$column];
	}
    }
    close($FILEHANDLE);
    print "ready.\n";
    return @columncontents;
}

sub convert_time_to_sec{
    # convert time elapsed to seconds
    # usage: $sec=convert_time_to_sec($day1,$day2,$hour1,$hour2,$min2,$min1,$sec1,$sec2)
    my $secs=0;
    if ($sec2<$sec1) {
	$sec2+=60;
	$min2--;
    }
    $secs+=($sec2-$sec1);
    if ($min2<$min1) {
	$min2+=60;
	$hour2--;
    }
    $secs+=(60*($min2-$min1));
    if ($hour2<$hour1) {
	$hour2+=24;
	$day2--;
    }
    $secs+=(3600*($hour2-$hour1));
    $secs+=86400*($day2-$day1);
}

1;   #it took me about one hour to realise i had to put this guy here!
