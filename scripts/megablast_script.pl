use warnings;
use strict;

# This script let's us run megablast successively to compare each .fasta result file with the others.

# It should :
# * Take arguments on command line
# * Inject them into the template database and megablast commands variables
# * Execute said template variables using each .fasta file in a given directory

# It should do this using each .fasta as a database, and megablast the others against it.

# AKA., 2,3,4 vs DB1;  1,3,4 vs DB2;  1,2,4 vs DB3 ...

if (scalar @ARGV != 2){die "Program usage:\nmegablast_script.pl <Input_fasta_directory> <Results_directory>\n"}
my $input_dir=$ARGV[0];
my $out_dir=$ARGV[1];



my @liste_fasta = glob( $input_dir . '/*.fasta' );


my $db_subject_name;
my $db_subject_path;


# for each .fasta :
## Take 1 .fasta to serve as the DB,
## Run it against the database of each other .fasta in the folder,
## save each result to blastRes/megablast_Test/db_subject.csv

# take an element to serve as the DB
foreach my $db (@liste_fasta) {

  #take another to serve as the subject
  foreach my $subject (@liste_fasta){
    # if two elements identical do nothing
    # else check if el1_el2.csv exists in output dir:
    if ($db ne $subject){

      # extract the name of the db from  "directory/name.fasta" and store it in #db_id
      $db =~ /.*[\/](.*)[.].*/;
      my $db_id = $1;


      # extract the name of the subject from  "directory/name.fasta" and store it in #subject_id
      $subject =~ /.*[\/](.*)[.].*/;
      my $subject_id = $1;


      # create two variables. one containing the full to give to the resulting csv file, the other containing the path to the resulting file.

      # EXISTS:
      ## do nothing
      # DOESN'T EXIST:
      ## Create the new file name:
      ## Run megablast with el1 as DB and el2 as subject. name output el1_el2.csv


      $db_subject_name =  $db_id . "_vs_" . $subject_id . ".csv";
      ##$db_subject_path = "blastRes/megablast_Test/" . $db_subject_name;

      system("megablast -d $db -i $subject  -o $out_dir/$db_subject_name -q 5 -r -4 -D 3 -a 4 -W 32 -m 8");




    }
  }
}
