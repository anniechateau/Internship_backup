#!/usr/bin/perl
use strict;
use DBI;
use List::Util qw(any);

# This script parses the blast results and created a results matrix to create
# a

## Only executes the script if it receives the 3 necessary arguments:
## The matrice file in which to store the comparison matrice.
#  (filename to write matrice to. should be created by this program / rewritten if already existing)
## the output directory in which to place any files created.
#  (optional so as to store .dist files only when asked ?)
## The list directory where all the blast results are stored.
if (scalar @ARGV != 3){die "Program usage:\nparsing.pl <Matrice_file> <Output_directory> <blast_results_directory>\n"}
my $mat_file=$ARGV[0];
my $out_dir=$ARGV[1];
my $dir_csv=$ARGV[2];

## temporary value. to replace by function that determines the actual length of the original data, and divides by 2 (or user input).
## since we're using a fasta file.
my $nbre = 100;

# Script pour rentrer les résultats de la comparaison dans la matrice.
# liste des fichiers csv contenant les resultats des blasts.
my @liste_csv = glob( $dir_csv . '/*.csv' );
foreach(@liste_csv){

	# 1. read each blast output result file, summing scores for each one.

	# for each CSV file containing out blast results:
	print "\non traite $_\n";
	# for every CSV, extract the ID's from the name, assuming it is named:
	# ID1-vs-ID2.csv
	my $fic=$_;
	$fic=~/.*\/([^\/]*)-vs-([^\/]*)\.csv/;
	my $id1=$1;
	my $id2=$2;
	# read the result and get the score for each line.
	open (F1,"<$fic") || die "\npblm ouverture $fic\n";
	my %vu_f;
	my $score_id1_id2=0;
	my $n_b=0;
	my $check_blast=0;
	while (<F1>){
		my @ligne=split(",",$_);
		if ($#ligne!=14 & $#ligne!=0){$check_blast=1;}
		if((!defined($vu_f{$ligne[0]})) && ($ligne[0] ne $ligne[1])){
			$vu_f{$ligne[0]}=1;
			$score_id1_id2+=$ligne[11];
			#print "ligne 11 $ligne[11]\n";
			$n_b++;
		}
	}
	close F1;


## This code exists to relaunch every Blast comparison in case there's some error in it
## (in this case, error being one of the line's not having 14 values in it).
##
## I think it'd be best to create a pipeline, where one of the amont actions launches the BLASTS and stores the results.
	if ($check_blast==1){
		die "\n###!!!### pblm avec blast $id1 vs $id2 -> on relance le tout\n";
	}
	else{
		if ($n_b==0){
			# die "pblm entre $id1 et $id2, pas de résultat de blast\n";
			$score_id1_id2="0.0000";
		}
		$score_id1_id2/=$nbre;
		print "$id1 vs $id2 : $n_b\t$score_id1_id2\n";
		## The global score of one of the blast files is considered to be:
		## The summed score of every alignment, divided by the total number of sequences/alignments.


		## 2. Save scores for each file in .dist files (distance files) with the following syntax:

		# id1_id2.dist contains (separated by tabs)
		# "id1	id2		id1-id2-score"


		my $out_file_id1_id2=$out_dir."/".$id1."_vs_".$id2.".dist";
		my $out_temp = $out_dir."/".$id1."_vs_".$id2.".dist_temp";
		if (-e $out_file_id1_id2){
			open(DIST,"<$out_file_id1_id2") || die ("pblm opening file $out_file_id1_id2\n");
			open(TEMP,">$out_temp") || die ("pblm opening file $out_temp\n");
			while(<DIST>){
				chomp($_);
				my @csv=split(",",$_);
				if ($csv[0]==$id2){
					print "Already a result for $id1 vs $id2\n";
				}
				else{
					print TEMP "$_\n";
				}
			}
			print TEMP "$id1\t$id2\t$score_id1_id2\n";
			close DIST;
			close TEMP;
			`mv $out_temp $out_file_id1_id2`;
		}
		else{
			print "no file $out_file_id1_id2 - first comparison -> we create the file\n";
	                open(DIST,">$out_file_id1_id2") || die ("pblm opening file $out_file_id1_id2\n");
			print DIST "$id1\t$id2\t$score_id1_id2\n";
			close DIST;

		}
	}
}



 #  3. Create Matrix csv file:
 #if (any { $_ eq $id1 } @matrice){
 #	print @matrice;
 #	my ;
 #	push(@matrice, \@id1);
 #}


##for each dist file:
# add to @list_ids if not already present.

my @list_ids;
# get all .dist files in output directory.
my @list_dists = glob( $out_dir . '/*.dist' );
foreach(@list_dists){
	# take current filename and test it against the RegEx
	my $query=$_;
	$query=~/.*\/([^\/]*)_[^\/]*\.dist/;
	# extract id1 from the RegEx.
	my $newID = $1;

	# inefficient, need to check if there's a better, faster, way to see if an item is already in an array.
	if (any { $_ eq $newID } @list_ids){
		# pass, as we don't want several copies of the ID in our
	}else {
		push(@list_ids,$newID);
	}
}



my $line = "";
my $currentID;
my $currentDistFile;
my $score;
my $query;
my $otherID;


#open matrice.csv
open(MAT,">$mat_file") || die ("pblm opening file $mat_file\n");


# write " ,id1,id2,id3,[...]" with all id's in array, as first line of matrice csv file.
print MAT " ";
foreach(@list_ids){
	print MAT ",".$_;
}
print MAT "\n";




for my $i (0 .. scalar @list_ids -1){

	$currentID = $list_ids[$i];
	$line = $currentID . ",";
	print "test Line should contain current ID -> $line end line";
	print "\n\ncurrent main ID: $currentID\n";

	# iterate over ID' in that array again.
	for my $j (0 .. scalar @list_ids -1){
		# if ID outer loop == ID' inner loop:
		$otherID = $list_ids[$j];
		print "\n\ncurrent secondary ID: $otherID\n";

		# add comma separator if this isn't the first item in the list.


		if($i == $j){
			# if we're at the intersection point (AKA., A vs A, B vs B...)
			# then set the value to 0. (or should it be set to a very high score ?)
			$line = $line . "0";
			print "\n\tSame ID's --> write 0\n";
			}
		elsif($i != $j) {
			# Get the current distance file for ID1 ($currentID) vs ID2 ($otherID):
			# open the file in Read mode.
			$currentDistFile = $out_dir."/".$currentID."_".$otherID.".dist";
			print "\n\tdifferent ID's --> read $currentDistFile\n";
			open(DISTFILE,"<$currentDistFile") || die ("pblm opening file $currentDistFile for querying\n");
			while(<DISTFILE>){
				$query = $_;
				chomp $query;
				print "\n\t\tQuery is $query\n";
				 #example subject: BM-CLC	BM-MegaHit	(18.418)   <--- we want the score.

				$query=~/\w*-\w*\s\w*-\w*\s(\d+.?\d*)/;
				$line = $line . "$1";
				print "\n\t\tRecognized pattern -> $1\n";
			}
			close(DISTFILE);

		}





		# write the generated line into the matrice file.
		if ($j != $#list_ids){
			$line = $line . ",";
		}
	}
	print "\nFinal Line is: $line\n";
	print MAT $line . "\n";

}
			#  line += "0,"
		# else:
			# read distance score file ID_ID'.dist
			# extract score
			# line += "score,"   # inner loop ID', not outer loop ID
	# write $line to matrice.csv

#close matrice.csv
