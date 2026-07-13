#!/usr/local/perl

use strict;
use FileHandle;

return 1;

# this procedure takes the "results" output of score.c (when MESSAGE is defined)
# and sorts the different categories of splice sites after gene_id and then after
# location of splice site in the gene
sub order {
    my ($input,$output)=@_;

    open(Inp,$input) or die "Error 81: coudn't open $input for reading!\n";
    open(Out,">$output") or die "Error 82: couldn't open $output for writing!\n";

    my $status;
    my $o="temp.out";

    while (<Inp>) {

	if($_ eq "Scores for false donors\n") {
	    O->autoflush(1);
	    close(O);
	    $status=system("sort -k 1,2 temp.out>temp.sort.out");
	    die "ERROR 83: Could not sort temporary file when filtering: $!\n" unless $status == 0;
	    open(F,"temp.sort.out") or die "Error 84: couldn't open temporary file temp.sort.out for reading: $!\n";
	    while(<F>) { print Out $_;}
	    close(F);
	    print  Out "Scores for false donors\n";
	    open(O,">$o") or die "Error 85: couldn't open temporary file $o for writing: $!\n";
	}
	elsif($_ eq "Scores for false acceptors\n"){
	    O->autoflush(1);
	    close(O);
	    $status=system("sort -k 1,2 temp.out>temp.sort.out");
	    die "ERROR 86: Could not sort temporary file when filtering: $!\n" unless $status == 0;
	    open(F,"temp.sort.out") or die "Error 87: couldn't open temporary file temp.sort.out for reading: $!\n";
	    while(<F>) { print Out $_;}
	    print  Out "Scores for false acceptors\n";
	    close(F);
	    open(O,">$o") or die "Error 88: couldn't open temporary file $o for writing: $!\n";
	}
	elsif($_ eq "Scores for true acceptors\n"){
	    print Out "Scores for true acceptors\n";
	    open(O,">$o") or die "Error 89: couldn't open temporary file $o for writing: $!\n";
	}
	elsif($_ eq "Scores for true donors\n"){
	    O->autoflush(1);
	    close(O);
	    $status=system("sort -k 1,2 temp.out>temp.sort.out");
	    die "ERROR 90: Could not sort temporary file when filtering: $!\n" unless $status == 0;
	    open(F,"temp.sort.out") or die "Error 91: couldn't open temporary file temp.sort.out for reading: $!\n";
	    while(<F>) { print Out $_;}
	    print Out "Scores for true donors\n";
	    close(F);
	    open(O,">$o") or die "Error 92: couldn't open temporary file $o for writing: $!\n";
	}
	else {
	    my @a=split(/ +/);
	    print O $a[0]," ";
	    printf O "%010d ",$a[1];
	    print O "@a[2..$#a]";
	}
    }

    O->autoflush(1);
    close(O);
    $status=system("sort -k 1,2 temp.out>temp.sort.out");
    die "ERROR 93: Could not sort temporary file when filtering: $!\n" unless $status == 0;
    open(F,"temp.sort.out") or die "Error 94: couldn't open temporary file temp.sort.out for reading: $!\n";
    while(<F>) { print Out $_;}
    close(F);

    Out->autoflush(1);
    close(Out);
    close(Inp);

    system("rm temp.out temp.sort.out");
}


# this procedure takes the ordered results output of score program when
# MESSAGE is defined and keeps only the best scores within a given
# distance (usually 30bp or 60 bp)
sub eliminsites {
    my ($input,$output,$filtacc,$filtdon)=@_;

    my $dist=$filtacc;

    open(F,$input) or die "Error 95: coudn't open $input for reading: $!\n";
    open(Out,">$output") or die "Error 96: couldn't open $output for writing: $!\n";

    my $kt=0;
    my $kf=0;
    my $ta=0;
    my $td=0;
    my $fa=0;
    my $fd=0;

    my (@scoret,@scoret1,@scoret2,@scoret3,@genet,@pozt);
    my (@scoref,@scoref1,@scoref2,@scoref3,@genef,@pozf,%isfa,%isfd);

    while(<F>) {
    
	chomp;

	my @a=split;

	if($ta && $#a == 5) {
	    $scoret[$kt]=$a[5];
	    $scoret1[$kt]=$a[2];
	    $scoret2[$kt]=$a[3];
	    $scoret3[$kt]=$a[4];
	    $genet[$kt]=$a[0];
	    $pozt[$kt++]=$a[1];
	}

	if($fa && $#a == 5) {
	    $scoref[$kf]=$a[5];
	    $scoref1[$kf]=$a[2];
	    $scoref2[$kf]=$a[3];
	    $scoref3[$kf]=$a[4];
	    $genef[$kf]=$a[0];
	    $pozf[$kf++]=$a[1];
	    $isfa{$a[0]}++;
	}

	if($td && $#a == 5) {
	    $scoret[$kt]=$a[5];
	    $scoret1[$kt]=$a[2];
	    $scoret2[$kt]=$a[3];
	    $scoret3[$kt]=$a[4];
	    $genet[$kt]=$a[0];
	    $pozt[$kt++]=$a[1];
	}

	if($fd && $#a == 5) {
	    $scoref[$kf]=$a[5];
	    $scoref1[$kf]=$a[2];
	    $scoref2[$kf]=$a[3];
	    $scoref3[$kf]=$a[4];
	    $genef[$kf]=$a[0];
	    $pozf[$kf++]=$a[1];
	    $isfd{$a[0]}++;
	}

	if($_ eq "Scores for true acceptors") {
	    $ta=1;
	}

	if($_ eq "Scores for false acceptors") {
	    $ta=0;
	    $fa=1;
	}

	if($_ eq "Scores for true donors") {

	    my ($gene,$e1,$b1,$b2,$e2,$e2,$u,$last,$i,$i1,$i2,$val,@rest,@score,@score1,@score2,@score3,@ind,@elim);

	    # see the missing sites
	    print Out "Acceptors:\n";
	    $b1=0;
	    $b2=0;
	    $u=0;

	    while($b1<$kt) {
		$gene=$genet[$b1];
		$e1=$b1;
		while(($genet[$e1] eq $gene)&&($e1<$kt)) { $e1++;}
		if($isfa{$gene}) {
		    $e2=$b2;
		    while(($genef[$e2] eq $gene)&&($e2<$kf)) { $e2++;}

		    if($e2==$b2) {
			$gene=$genef[$b2];
			$e2=$b2;
			while(($genef[$e2] eq $gene)&&($e2<$kf)) { $e2++;}
			$i2=$b2;
			$last=$u;		
			for($i=$last;$i<$last+$e2-$b2;$i++) {
			    $rest[$i]=$pozf[$i2];
			    $score[$i]=$scoref[$i2];
			
			    $score1[$i]=$scoref1[$i2];
			    $score2[$i]=$scoref2[$i2];
			    $score3[$i]=$scoref3[$i2];

			    $i2++;
			    $ind[$i]=0;
			    $elim[$i]=0;
			
			    while($rest[$i]>$rest[$u]+$dist) { $u++;}
			    $val=$score[$i];

			    for(my $j=$u;$j<$i+1;$j++) {
				if(($score[$j]<$val) && (!$elim[$j])) { $elim[$j]=1;}
				else { if($score[$j]>$val) {$val=$score[$j];}}
			    }
			}

			$b2=$e2;
			$u=$i;
			
		    }	    
		    else {
		    
			$i1=$b1;
			$i2=$b2;
			$last=$u;
		
			for($i=$last;$i<$last+$e1-$b1+$e2-$b2;$i++) {
			
			    if(($i1<$e1 && $pozt[$i1]<$pozf[$i2]) || ($i2==$e2 && $i1<$e1) ) {
				$rest[$i]=$pozt[$i1];
				$score[$i]=$scoret[$i1];

				$score1[$i]=$scoret1[$i1];
				$score2[$i]=$scoret2[$i1];
				$score3[$i]=$scoret3[$i1];

				$ind[$i]=1;
				$elim[$i]=0;
				$i1++;
			    
			    }
			    else {
				$rest[$i]=$pozf[$i2];
				$score[$i]=$scoref[$i2];
				
				$score1[$i]=$scoref1[$i2];
				$score2[$i]=$scoref2[$i2];
				$score3[$i]=$scoref3[$i2];
				
				$i2++;
				$ind[$i]=0;
				$elim[$i]=0;			
			    
			    }

			    while($rest[$i]>$rest[$u]+$dist) { $u++;}
			    $val=$score[$i];

			    for(my $j=$u;$j<$i+1;$j++) {
				if($score[$j]<$val && !$elim[$j]) { $elim[$j]=1;}
				else { if($score[$j]>$val) {$val=$score[$j];}}
			    }
			}

			$b1=$e1;
			$b2=$e2;
			$u=$i;
		    }
		}
		else { 
		    $b1=$e1;
		}
	    }

	    # print the results
	    for(my $i=0;$i<$kt+$kf-1;$i++) {

		print Out $score[$i]," ",$ind[$i]," ",$rest[$i]," ",$elim[$i];
		if($ind[$i] && $elim[$i] ) {print Out " *";}
		print Out "\n";
	    }

	    $fa=0;
	    $td=1;
	    $kt=0;
	    $kf=0;
	}

	if($_ eq "Scores for false donors") {
	    $td=0;
	    $fd=1;
	}

    }

    $dist=$filtdon;

    # see the missing sites
    print Out "Donors:\n";

    my ($gene,$e1,$b1,$b2,$e2,$e2,$u,$last,$i,$i1,$i2,$val,@rest,@score,@score1,@score2,@score3,@ind,@elim);

    $b1=0;
    $b2=0;
    $u=0;

    while($b1<$kt) {

	$gene=$genet[$b1];
	$e1=$b1;
	while(($genet[$e1] eq $gene)&&($e1<$kt)) { $e1++;}
	if($isfd{$gene}) {
	    $e2=$b2;
	    while(($genef[$e2] eq $gene)&&($e2<$kf)) { $e2++;}

	    if($e2==$b2) {

		$gene=$genef[$b2];
		$e2=$b2;
		while(($genef[$e2] eq $gene)&&($e2<$kf)) { $e2++;}
		$i2=$b2;
		$last=$u;		
		for($i=$last;$i<$last+$e2-$b2;$i++) {
		    $rest[$i]=$pozf[$i2];
		    $score[$i]=$scoref[$i2];

		    $score1[$i]=$scoref1[$i2];
		    $score2[$i]=$scoref2[$i2];
		    $score3[$i]=$scoref3[$i2];
		    
		    $i2++;
		    $ind[$i]=0;
		    $elim[$i]=0;
	    
		    while($rest[$i]>$rest[$u]+$dist) { $u++;}
		    $val=$score[$i];

		    for(my $j=$u;$j<$i+1;$j++) {
			if(($score[$j]<$val) && (!$elim[$j])) { $elim[$j]=1;}
			else { if($score[$j]>$val) {$val=$score[$j];}}
		    }
		}

		$b2=$e2;
		$u=$i;
	
	    }	    
	    else {

		$i1=$b1;
		$i2=$b2;
		$last=$u;
		
		for($i=$last;$i<$last+$e1-$b1+$e2-$b2;$i++) {
	    
		    if(($i1<$e1 && $pozt[$i1]<$pozf[$i2]) || ($i2==$e2 && $i1<$e1) ) {
			$rest[$i]=$pozt[$i1];
			$score[$i]=$scoret[$i1];

			$score1[$i]=$scoret1[$i1];
			$score2[$i]=$scoret2[$i1];
			$score3[$i]=$scoret3[$i1];

			$ind[$i]=1;
			$elim[$i]=0;
			$i1++;

		    }
		    else {
			$rest[$i]=$pozf[$i2];
			$score[$i]=$scoref[$i2];
			
			$score1[$i]=$scoref1[$i2];
			$score2[$i]=$scoref2[$i2];
			$score3[$i]=$scoref3[$i2];

			$i2++;
			$ind[$i]=0;
			$elim[$i]=0;			
		    
		    }

		    while($rest[$i]>$rest[$u]+$dist) { $u++;}
		    $val=$score[$i];

		    for(my $j=$u;$j<$i+1;$j++) {
			if($score[$j]<$val && !$elim[$j]) { $elim[$j]=1;}
			else { if($score[$j]>$val) {$val=$score[$j];}}
		    }
		}

		$b1=$e1;
		$b2=$e2;
		$u=$i;
	    }
	}
	else {
	    $b1=$e1;
	}
    }

   # print the results
    for(my $i=0;$i<$kt+$kf-1;$i++) {

	print Out $score[$i]," ",$ind[$i]," ",$rest[$i]," ",$elim[$i];
	if($ind[$i] && $elim[$i] ) {print Out " *";}
	print Out "\n";
    }

    close(F);
    Out->autoflush(1);
    close(Out);
}


# this procedure takes the output of eliminsites and makes 4 
# files with the remained splice site scores
sub printres {
    
    my $res=$_[0];

    open(F,$res) or die "Error 97: coudn't open $res for reading: $!\n";

    my $o1="res.acc";
    my $o2="res.facc";
    my $o3="res.don";
    my $o4="res.fdon";

    my ($acc,$don);

    while(<F>)
    {
	chomp;

	my @a=split;

	if($_ eq "Acceptors:") { $acc=1; open(O1,">$o1"); open(O2,">$o2"); } else {
    
	    if($_ eq "Donors:") { 
		$don=1; $acc=0; 
		O1->autoflush(1);close(O1); 
		O2->autoflush(1);close(O2); 
		open(O3,">$o3");open(O4,">$o4");
	    } else {
	
		if($acc) {
		    if($a[3]) {
			if($a[1]) { print O1 "-99.00\n";}
			else { print O2 "-99.00\n";}
		    }
		    else { 
			if($a[1]) { print O1 $a[0],"\n";}
			else { print O2 $a[0],"\n";}
		    }
		}
		
		if($don) {
		    if($a[3]) {
			if($a[1]) { print O3 "-99.00\n";}
			else { print O4 "-99.00\n";}
		    }
		    else { 
			if($a[1]) { print O3 $a[0],"\n";}
			else { print O4 $a[0],"\n";}
		    }
		}
	    }
	}
    }

    close(F);

    O3->autoflush(1);close(O3); 
    O4->autoflush(1);close(O4); 
}

sub sortres{
    my ($input,$o,$message)=@_;

    open(O,">> $o") or die "Error 98: coudn't open $o for writing: $!\n";
    print O "$message\n";
    my $status=system("sort -n $input>temp.out");
    die "ERROR 99: Could not sort temporary file for printing!\n" unless $status == 0;
    open(F,"temp.out") or die "Error 100: coudn't open temporary file for reading!\n";
    while(<F>){
	print O $_;
    }
    close(F);
    
    O->autoflush(1);
    close(O);

    my $file="temp.out";
    unlink($file);
}
