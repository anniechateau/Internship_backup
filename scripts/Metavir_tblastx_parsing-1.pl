#!/usr/bin/perl
use strict;
use DBI;
# Script pour parser les résultats des blasts -> tourne depuis metagrid
# Constantes :
my $mat_file="metavir/Data/Virome_Comp/matrice_comp.tab";
my $out_dir="metavir/Virome_Comp/Data_storage/";
my $liste_dir="metavir/Virome_Comp/Lists/";
my $nbre=50000; # parce qu'on a des échantillons de 50 000 séquences

# Script pour rentrer les résultats de la comparaison dans la matrice.
# liste des fichiers tabs arrivés
my @liste_tab=<metavir/Virome_Comp/*.tab>;
foreach(@liste_tab){
	# pour chaque fichier tab
	print "on traite $_\n";
	my $fic=$_;
	$fic=~/.*\/([^\/]*)-vs-([^\/]*)\.tab/;
	my $id1=$1;
	my $id2=$2;
	# on parse le résultat
	# et on supprime le fichier tab
	open (F1,"<$fic") || die "pblm ouverture $fic\n";
	my %vu_f;
	my $score_id1_id2=0;
	my $n_b=0;
	my $check_blast=0;
	while (<F1>){
		my @ligne=split("\t",$_);
		if ($#ligne!=11){$check_blast=1;}
		if((!defined($vu_f{$ligne[0]})) && ($ligne[0] ne $ligne[1])){
			$vu_f{$ligne[0]}=1;
			$score_id1_id2+=$ligne[11];
			# print "ligne 11 $ligne[11]\n";
			$n_b++;
		}
	}
	close F1;
	if ($check_blast==1){
		print "###!!!### pblm avec blast $id1 vs $id2 -> on relance le tout\n";
		my $associated_tab_file="/storage/metavir/Virome_Comp/".$id2."-vs-".$id1.".tab";
#		my $file_sh_associee="/storage/metavir/scripts/Virome_comp/scripts_sh_temp/blast-".$id1."-".$id2.".sh";
#		my $file_sh_associee_cluster="/home/servers/metavir/scripts/Virome_comp/scripts_sh_temp/blast-".$id1."-".$id2.".sh";
#		if (-e $file_sh_associee){}
#		else{
#			$file_sh_associee="/storage/metavir/scripts/Virome_comp/scripts_sh_temp/blast-".$id2."-".$id1.".sh";
#			$file_sh_associee_cluster="/home/servers/metavir/scripts/Virome_comp/scripts_sh_temp/blast-".$id2."-".$id1.".sh";
#		}
		my $out="";
		if (-e $associated_tab_file){$out=`mv $associated_tab_file /storage/metavir/Virome_Comp/Pblmatic_files/`;}
		print "we move $associated_tab_file to the pblmatic files directory : $out\n";
		$out=`mv $fic /storage/metavir/Virome_Comp/Pblmatic_files/`;
		print "we move $fic to the pblmatic files directory : $out\n";
	}
	else{
		if ($n_b==0){
			# die "pblm entre $id1 et $id2, pas de résultat de blast\n";
			$score_id1_id2="0.0000";
		}
		$score_id1_id2/=$nbre;
		print "$id1 vs $id2 : $n_b\t$score_id1_id2\n";
		# on remplit les fichiers dist en mettant a jour la matrice
		my $out_file_id1=$out_dir.$id1.".dist";
		my $out_temp=$out_dir.$id1.".dist_temp";
		if (-e $out_file_id1){
			open(DIST,"<$out_file_id1") || die ("pblm opening file $out_file_id1\n");
			open(TEMP,">$out_temp") || die ("pblm opening file $out_temp\n");
			while(<DIST>){
				chomp($_);
				my @tab=split("\t",$_);
				if ($tab[0]==$id2){
					print "Already a result for $id2 vs $id1, we remove\n";
				}
				else{
					print TEMP "$_\n";	
				}
			}
			print TEMP "$id2\t$score_id1_id2\n";
			close DIST;
			close TEMP;
			`mv $out_temp $out_file_id1`;
		}
		else{
			print "no file $out_file_id1 - first comparison -> we create the file\n";
	                open(DIST,">$out_file_id1") || die ("pblm opening file $out_file_id1\n");
			print DIST "$id2\t$score_id1_id2\n";
			close DIST;

		}
	}
}
