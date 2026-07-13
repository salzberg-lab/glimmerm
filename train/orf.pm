#!/usr/local/perl

# this library takes the files seqs and exons.dat and extracts from them
# two output files:
# one with the orfs
# the other one with introns --- commented

use strict;
use FileHandle;

return 1;

sub genorf{ 
    my ($f,$g) = @_;

    my $f1="orfs";
#    my $f2="introns";

    my (%isexon, %isanum, %isintron, $end, $start, $start1, $end1, $anum);

    open(F,$f) or die "ERROR 131: Couldn't open file $f for reading!\n";

    while(<F>){
	chomp;
	($anum,$start,$end) = split;

	if($isanum{$anum}) {
	    $isexon{$anum}.="$start $end ";
	    $start1=$start-1;
	    $isintron{$anum}.="$end1 $start1 ";	
	}

	else {
	    if($anum ne "") {
		$isanum{$anum}=1;
		$isexon{$anum}="$start $end ";
	    }
	}

	$end1=$end+1;
	
    }
    close(F);
	

    open(G,$g) or die "ERROR 132: Couldn't open file $g for reading!\n";

    open(F1,">$f1") or die "ERROR 133: Couldn't open file $f1 for writing!\n";
#    open(F2,">$f2");

    my $seq;

    while(<G>){
	chomp;
	($anum,$seq)=split;
    
	my @e=split(/\s+/,$isexon{$anum});
	my @i=split(/\s+/,$isintron{$anum});

	my $orf="";

	for(my $k=0;$k<=$#e;$k+=2){
	    $orf.=substr($seq,$e[$k]-1,$e[$k+1]-$e[$k]+1);
	}

	# my comment begin
	# if you want the whole coding gene use:
	# print F1 "$anum $orf\n";

	# otherwise: skip start condon and last codon:
	
	my $l=length($orf);
	my $neworf=substr($orf,3,$l-6);

	print F1 "$anum $neworf\n";

	# my comment end

#	for(my $k=0;$k<$#i;$k+=2){
#	    my $intron=substr($seq,$i[$k]-1,$i[$k+1]-$i[$k]+1);
#	    print F2 "$anum $intron\n";
#	}
    }

    close(G);
    
    F1->autoflush(1);
    close(F1);

    STDOUT->autoflush;
}

