#!/usr/bin/perl
# polarize_reads.pl script. Custom script used in the curation pipeline for PacBio reads as descibed in Jamy et al. 2019
# By Mahwash

use strict;
use warnings;

die "\nReads in barrnap outfile and determines which reads need to be reverse-complemented. Outputs fasta file with polarized reads.\n\nUsage: polarize_reads.pl <barrnap file> <fasta file> <output>\n\n" unless @ARGV == 3;

my ($barrnap, $fasta, $output) = @ARGV;

my %fasta = ();
my $header = "";
my $concatenated_seq = "";
my %read_dir = ();
my @columns = ();
my $dir = "";

# open barrnap file and read line by line. Create hash which stores direction for each seq id. Barrnap file contains several lines for each seq id, corresponding to 18S, 28S and sometimes 5.8S. To make things easier, only look at lines containing 18S.
open (my $in_barrnap, "<$barrnap") or die "error opening $barrnap for reading";
while (my $line = <$in_barrnap>) {
    if ($line =~ /18S_rRNA/) {
            @columns = split("\t", $line);
            ($header, $dir) = ($columns[0], $columns[6]);
            $read_dir{$header} = $dir;
    }
}
        
close $in_barrnap;


# Build a hash where the keys (headers in fasta file) point at the sequence
open(my $in_fasta, "<$fasta") or die "error opening $fasta for reading";

while (my $line = <$in_fasta>) {
    if ($line =~ /^>(\s*\S+).*/) {
        $header = $1;
        $concatenated_seq = "";
    }
    else {
        chomp($line);
        $concatenated_seq = $concatenated_seq . $line;
        $fasta{$header} = $concatenated_seq;
    }
}

close $in_fasta;


# To outfile, print 'positive' reads as they are and reverse complement 'negative' reads.
open (my $out, ">$output") or die "error opening $output for writing";

foreach my $header (keys %fasta) {
    if ($read_dir{$header} eq "+" ) {
        print $out ">$header\n$fasta{$header}\n";
    } elsif ($read_dir{$header} eq "-") {
    	print $out ">$header\n"; 
        my $revcomp = reverse($fasta{$header});
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;        
        print $out "$revcomp\n";
    }       
}               

close $out;
