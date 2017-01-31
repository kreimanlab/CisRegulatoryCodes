#!/usr/bin/perl
# retrieve nm/ug/ll from any of these variables

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";

$n_args=$#ARGV+1;
if ($n_args < 4) {
    print "ll2nm2ug.pl <search_filename> <search_col> <species> <searchtype>\n";
    print "output is stored in query.txt (hits) and log.txt (miss+multiple matches)\n";
    print "searchtype is ll / nm / ug\n";
    gk();
    exit;
}

$search_filename=$ARGV[0];
$search_col=$ARGV[1];
$species=$ARGV[2];
$search_col--;
$searchtype=$ARGV[3];
#if ($searchtype eq "ll") {
#    $searchtype=1;
#} else {
#    if ($searchtype eq "nm") {
#	$searchtype=2;
#    } else {
#	if ($searchtype eq "ug") {
#	    $searchtype=3;
#	}
#    }
#}
print "search_filename=$search_filename\n";
print "search_col=$search_col\n";
print "species=$species\n";
print "searchtype=$searchtype\n";

if ($species eq "hs") {
    @homol_seq_pairs_ll=read_column("${DATA_DIR}/db/locus_link/homol_seq_pairs_hsmm.txt",2,1,1);
} else {
    @homol_seq_pairs_ll=read_column("${DATA_DIR}/db/locus_link/homol_seq_pairs_hsmm.txt",7,1,1);
}    
@homol_seq_pairs=read_whole_file_v2("${DATA_DIR}/db/locus_link/homol_seq_pairs_hsmm.txt",1,0);
if ($species eq "hs") {
    @homol_seq_pairs_ll2=read_column("${DATA_DIR}/db/homology/hsmm.txt",1,1,1);
} else {
    @homol_seq_pairs_ll2=read_column("${DATA_DIR}/db/homology/hsmm.txt",3,1,1);
}
@evidence=read_column("${DATA_DIR}/db/homology/hsmm.txt",5,1,0);
@percentages=read_column("${DATA_DIR}/db/homology/hsmm.txt",6,1,0);

@homol_seq_pairs2=read_whole_file_v2("${DATA_DIR}/db/homology/hsmm.txt",1,0);
@loc2ref=read_whole_file("${DATA_DIR}/db/locus_link/loc2ref");
@loc2UG=read_whole_file("${DATA_DIR}/db/locus_link/loc2UG");
#if ($species eq "mm") {
#    $map_filename="LL3_030128.mm.txt";
#} else {
#    $map_filename="LL3_030128.hs.txt";
#}
#@llnmug=read_whole_file($map_filename); # ll nm ug chr
#@querymap=read_column($map_filename,$searchtype,0,1);

open (hits_file,">query.txt") || die "could not open query.txt for output";
open (log_file,">log.txt") || die "could not open log.txt for output";
open (input_file,$search_filename) || die "could not open $search_filename for reading";

$i=0;
while ($record = <input_file>) {
    chomp $record;
    if ($record !~ /^>/) {
	@splitarray=split /\t/,$record;
	$curr_query=$splitarray[$search_col];

	$i++;
	print "i=$i\tquery=$curr_query...\n";
	$ll=-1;$nm=-1;$ug=-1;
	# try to get the locuslink
	if ($searchtype ne "ll") {
	    if ($searchtype eq "nm") {
		$nm=$curr_query;
		@searchstuff=grep /\b${curr_query}\b/,@loc2ref;
		if (!@searchstuff) {
		    print log_file "Error! i=$i\tquery=$curr_query\ti could not find the ll in loc2ref\n";
		} else {
		    $n_matches=$#searchstuff+1;
		    if ($n_matches>1) {
			print log_file "Warning! i=$i\tquery=$curr_query\tfound $n_matches matches in loc2ref\n";
		    }
		    $ll=get_column($searchstuff[0],0);
		}
		$ug=get_ug_from_ll($ll);
	    }
	    if ($searchtype eq "ug") {
		$ug=$curr_query;
		@searchstuff=grep /\b${curr_query}\b/,@loc2UG;
		if (!@searchstuff) {
		    print log_file "Error! i=$i\tquery=$curr_query\ti could not find the ll in loc2ug\n";
		} else {
		    $n_matches=$#searchstuff+1;
		    if ($n_matches>1) {
			print log_file "Warning! i=$i\tquery=$curr_query\tfound $n_matches matches in loc2ug\n";
		    }
		    $ll=get_column($searchstuff[0],0);
		}
		$nm=get_nm_from_ll($ll);
	    }
	} else {
	    $ll=$curr_query;
	    $ug=get_ug_from_ll($ll);
	    $nm=get_nm_from_ll($ll);
	}

	# try to get the homologous gene
	$ll2=-1;$nm2=-1;$ug2=-1;
	if ($ll ne "-1") {
	    @searchstuff=grep /\b${ll}$/,@homol_seq_pairs_ll;
	    if (@searchstuff) {
		$n_withhomologous++;
		$n_matches=$#searchstuff+1;
		if ($n_matches>1) {
		    print log_file "Warning i=$i i found $n_matches matches for ll=$ll in homol_seq_pairs_ll\n";
		}
		$ind=get_column($searchstuff[0],0);
		$line=$homol_seq_pairs[$ind];
		@splitarray=split /\t/,$line;
		if ($species eq "hs") {
		    $ll2=$splitarray[6];
		    $nm2=$splitarray[8];
		} else {
		    $ll2=$splitarray[1];
		    $nm2=$splitarray[3];
		}
		$ug2=get_ug_from_ll($ll2);
	    } else {
		print log_file "Warning! i=$i i could not find ll=$ll in homol_seq_pairs_ll\n";
		@searchstuff=grep /\b${ll}$/,@homol_seq_pairs_ll2;
		if (@searchstuff) {
		    $n_withhomologous++;
		    $n_matches=$#searchstuff+1;
		    if ($n_matches>1) {
			print log_file "Warning i=$i i found $n_matches matches for ll=$ll in homol_seq_pairs_ll2\n";
			($output_index,$best_evidence,$best_percentage)=chooseone(@searchstuff);
			print log_file "\tchose index=$output_index (best_evidence=$best_evidence, best_percentage=$best_percentage)\n";
			#print "searchstuff=@searchstuff\n";
			#print "output_index=$output_index\n";
			#print "best_evidence=$best_evidence\n";
			#print "best_percentage=$best_percentage\n";
			#exit;
		    } else {
			$output_index=0;
		    }
		    #$ind=get_column($searchstuff[0],0);
		    $ind=get_column($searchstuff[$output_index],0);
		    $line=$homol_seq_pairs2[$ind];
		    @splitarray=split /\t/,$line;
		    if ($species eq "hs") {
			$ll2=$splitarray[2];
			$ll2=~s /ll\.//;
			$ll2=~s /LL\.//;
			$nm2=$splitarray[3];
			if ($nm2 !~ /NM/) {
			    $nm2=get_nm_from_ll($ll2);
			}
		    } else {
			$ll2=$splitarray[0];
			$ll2=~s/ll\.//;
			$ll2=~s/LL\.//;
			$nm2=$splitarray[1];
			if ($nm2 !~ /NM/) {
			    $nm2=get_nm_from_ll($ll2);
			}
		    }
		    $ug2=get_ug_from_ll($ll2);
		} else {
		    print log_file "Warning! i=$i i could not find ll=$ll in homol_seq_pairs_ll2\n";
		}
	    }
	    $n_hits++;
	} else {
	    $n_miss++;
	}
	
	if ($ll ne "-1") {
	    $ll="ll.${ll}";
	}
	if ($ll2 ne "-1") {
	    $ll2="ll.${ll2}";
	}

	$nm=remove_version($nm);
	$nm2=remove_version($nm2);
	print hits_file "$i\t$ll\t$nm\t$ug\t$ll2\t$nm2\t$ug2\n";
    }
}
close (input_file);
close (hits_file);

print "n_processed=$i\n";
print "n_withhomologous=$n_withhomologous\n";
print "n_hits=$n_hits\n";
print "n_miss=$n_miss\n";
print log_file "n_processed=$i\n";
print log_file "n_withhomologous=$n_withhomologous\n";
print log_file "n_hits=$n_hits\n";
print log_file "n_miss=$n_miss\n";

close (log_file);

# end main

sub get_ug_from_ll
{
    # get ug from ll
    # $ug=get_ug_from_ll($ll);

    my @searchstuff;
    my $curr_ll=$_[0];
    my $ug;
    my $n_matches;
    my $curr_ug=-1;

    if (length($curr_ll)>0) {
	if ($curr_ll ne "-1") {
	    @searchstuff=grep /^${curr_ll}\b/,@loc2UG;
	    if (!@searchstuff) {
		print log_file "Warning! i could not found UG for curr_ll=$curr_ll\n";
	    } else {
		$n_matches=$#searchstuff+1;
		if ($n_matches>1) {
		    print log_file "Warning! i found $n_matches matches for curr_ll=$curr_ll in loc2UG\n";
		}
		$curr_ug=get_column($searchstuff[0],1);
	    }
	}
    } else {
	print "ll=$curr_ll\tlength<=0\n";
	exit;
    }
    return $curr_ug;
}
sub get_nm_from_ll
{
    # get nm from ll
    # $nm=get_nm_from_ll($ll);

    my @searchstuff;
    my $curr_ll=$_[0];
    my $nm;
    my $n_matches;
    my $curr_nm=-1;

    if ($curr_ll ne "-1") {
	@searchstuff=grep /^${curr_ll}\b/,@loc2ref;
	if (!@searchstuff) {
	    print log_file "Warning! i could not found nm for curr_ll=$curr_ll\n";
	} else {
	    $n_matches=$#searchstuff+1;
	    if ($n_matches>1) {
		print log_file "Warning! i found $n_matches matches for curr_ll=$curr_ll in loc2ref\n";
	    }
	    $curr_nm=get_column($searchstuff[0],1);
	}
    }

    return $curr_nm;
}

sub remove_version{
    # given an nm entry with a version number, remove the version
    my $old_nm=$_[0];
    my $new_nm;
    my $i;

    $i=index($old_nm,".");
    if ($i>=0) {
	$new_nm=substr($old_nm,0,$i);
    } else {
	$new_nm=$old_nm;
    }

    return $new_nm;
}

sub chooseone{
    # this subroutine is called if there are multiple matches while searching in homology/hsmm.txt
    # a single preferred entry is chosen based on the following preference order:
    # column 5, c (curated) > B (homology among 3 species) > b (homology between 2 species)
    # column 6, for two entries with B, choose the one with the higher percentage
    # usage:
    # ($output_index,$best_evidence,$best_percentage)=chooseone(@entry);
    #
    # uses @evidence and @percentages

    my @entries=@_;
    my $output_index=0;
    my $best_percentage=0;
    my $best_evidence=1;        # c=3, B=2 and b=1
    my $current_percentage;
    my $current_evidence;
    my $n_entries=$#entries;
    my $current_entry;
    my $current_entry_line;
    my @split_entry;
    my $i;

    #print "entries=@entries\n";

    for ($i=0;$i<=$n_entries;$i++) {
	$current_entry=$entries[$i];
	$current_entry_line=get_line_number($current_entry);
	#print "current_entry_line=$current_entry_line\n";

	$current_evidence=$evidence[$current_entry_line];
	#print "current_evidence=$current_evidence\n";
	if ($current_evidence eq "c") {
	    $current_evidence=3;
	} else {
	    if ($current_evidence eq "B") {
		$current_evidence=2;
	    } else {
		$current_eveidence=1;
	    }
	}
	if ($current_evidence>$best_evidence) {
	    $output_index=$i;
	    $best_evidence=$current_evidence;
	    $best_percentage=$percentages[$current_entry_line];
	} else {
	    if ($current_evidence==$best_evidence) {
		$current_percentage=$percentage[$curr_entry_line];
		if ($current_percentage > $best_percentage) {
		    $output_index=$i;
		    $best_percentage=$current_percentage;
		}
	    }
	}
    }

    #print "output_index=$output_index\n";
    #exit;
    return ($output_index,$best_evidence,$best_percentage);
}  
