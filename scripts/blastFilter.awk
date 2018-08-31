#! /usr/bin/awk -f

BEGIN {
  FS="\t"
  complete=0;
  inclusions=0;
  total=0;
  # here we remove the previous files, so as not to add into them the new results.

  system("rm graph/complete.csv || true");
  system("rm graph/inclusions.csv || true");


  # The "|| true" part is a shortcut to ignore the error that would occur
  # if the file graph/complete.csv doesn't exist.

}

/^[^#]/ {
  total++;
  if ( $4/(substr($1,index($1,"sequence_length=")+16)) > 0.95 && $3 > 80 && $11 < 0.000000000000001) {
    if ((substr($1,index($1,"sequence_length=")+16) / (substr($2,index($2,"sequence_length=")+16))) < 1.5 && (substr($1,index($1,"sequence_length=")+16) / (substr($2,index($2,"sequence_length=")+16))) > 0.5) {
        print $0 >> "graph/complete.csv";
        complete++;
    } else {
        print $0 >> "graph/inclusions.csv";
        inclusions++;
    }
  }
}

END {
  print "Complete :", complete, "\n";
  print "Inclusions :", inclusions, "\n";
  print "Total :", total, "\n\n";
  print "Complete %:", (complete/total) * 100 , "%\n";
  print "Inclusions %:", (inclusions/total) * 100, "%\n";

}
