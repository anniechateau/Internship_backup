use strict;
use warnings;
use File::Copy;
my $originalFasta = <ARGV>;
print($originalFasta);
copy($originalFasta,"./../tagged_fasta/".$originalFasta) or die "copy failed $!";

BEGIN { $^I = ""; }
while (defined($_ = <ARGV>)) {
    s/>/>{FASTA}-/g;
}
continue {
    die "-p destination: $!\n" unless print $_;
}
