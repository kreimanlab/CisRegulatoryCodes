#!/usr/bin/perl
# list_cooc_motifs.pl
# create a list of the motifs in the cooc directory

$CODE_DIR=$ENV{'CODEDIR'};
$DATA_DIR=$ENV{'DATADIR'};
require "${CODE_DIR}/perl/lib/fileio_methods_v1.pl";

$n_args=$#ARGV+1;
if ($n_args < 2) {
    print "usage:\n";
    print "list_cooc_motifs.pl <exp_name> <exp_id> <order_constraint> <min_transcripts> <source>\n";
    exit;
}

$order_constraint=-1;
$min_transcripts=4;
$max_ul=5000;$max_el=5000;$max_il=5000;
$source="ncbi";
$init_max_d=1;$finit_max_d=4;$max_d_array[1]=25;$max_d_array[2]=50;$max_d_array[3]=100;$max_d_array[4]=200;
$init_max_n_motifs=2;$finit_max_n_motifs=4;
$threshold_column=38;

$curr_arg=0;$exp_name=$ARGV[$curr_arg];
$curr_arg++;
if ($n_args>$curr_arg) {
    $exp_id=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $order_constraint=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $min_transcripts=$ARGV[$curr_arg];
}
$curr_arg++;
if ($n_args>$curr_arg) {
    $source=$ARGV[$curr_arg];
}

print "exp_name=$exp_name\n";
print "exp_id=$exp_id\n";
$i=0;
for ($index_max_d=$init_max_d;$index_max_d<=$finit_max_d;$index_max_d++) {
    $max_d=$max_d_array[$index_max_d];
    for ($max_n_motifs=$init_max_n_motifs;$max_n_motifs<=$finit_max_n_motifs;$max_n_motifs++) {
	$filename_cooc_out="${DATA_DIR}/databh/${exp_name}/analysis/cooc/${exp_id}.${max_d}.${min_transcripts}.${max_n_motifs}.${order_constraint}.${threshold_column}.${max_ul}.${max_el}.${max_il}.cooc-${source}.out.txt";
	if (-e $filename_cooc_out) {
	    ($nt,$nc)=line_count($filename_cooc_out);$n_cooc=$nt-$nc;
	    if ($n_cooc>0) {
		open (input_file,$filename_cooc_out) || die "could not open $filename_cooc_out for reading";
		while ($record=<input_file>) {
		    chomp $record;
		    if (is_comment($record) ne "true") {
			@splitarray=split /\t/,$record;
			for ($j=0;$j<$max_n_motifs;$j++) {
			    $i++;
			    $allmotifs[$i]=$splitarray[$j];
			}
		    }
		}
		close (input_file);
	    }
	} else {
	    print "$filename_cooc_out does not exist\n";
	}
    }
}

open (output_file,">list_cooc_motifs.txt") || die "could not open list_cooc_motifs.txt for output";
if ($i>0) {
    @sortedlist = sort {$a <=> $b} @allmotifs;
    $previous_motif=-1;
    for ($j=1;$j<=$i;$j++) {
	if ($sortedlist[$j] ne $previous_motif) {
	    $previous_motif=$sortedlist[$j];
	    if (length($previous_motif) > 0 ) {
		print output_file "$previous_motif\n";
	    }
	}
    }
}
close (output_file);
