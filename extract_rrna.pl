#!/usr/bin/perl
# extract_rrna.pl. Custom script used in the curation pipeline for PacBio reads as descibed in Jamy et al. 2019
# By Mahwash

use strict;
use warnings;

die "\nTakes in fasta file with PacBio reads and extracts the 18S and 28S genes based on barrnap predictions as well as BLAST results in the case of 18S.\n\nUsage: extract_rrna.pl <barrnap file> <blast file> <fasta file> <18S output> <28S output>\n\n" unless @ARGV == 5;

my ($barrnap, $blast, $fasta, $SSUout, $LSUout) = @ARGV;

my %fasta = ();
my $header = "";
my $concatenated_seq = "";
my %blast_coords = ();
my @columns = ();
my $start = ();
my $end = ();
my $dir = ();
my %barrnap_18S = ();
my %barrnap_28S = ();
my $temp = "temp.txt";
my %extract = ();


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


# Build a hash where the keys (seq_id in blast file) point at an array containing the start and end coordinates
open (my $in_blast, "<$blast") or die "error opening $blast for reading";
while (my $line = <$in_blast>) {
    @columns = split("\t", $line);
    ($header, $start, $end) = ($columns[0], ($columns[6] - 1), $columns[7]);
    $blast_coords{$header} = [$start, $end];
}

close $in_blast;


# Build a hash where the keys (seq_id in barrnap file) point at an array containing the start and end coordinates
open (my $in_barrnap, "<$barrnap") or die "error opening $barrnap for reading";
while (my $line = <$in_barrnap>) {
    if ($line =~ /18S_rRNA/) {
        @columns = split("\t", $line);
        ($header, $start, $end, $dir) = ($columns[0], ($columns[3] - 1), $columns[4], $columns[6]);
        $barrnap_18S{$header} = [$start, $end, $dir];
    }elsif ($line =~ /28S_rRNA/) {
        @columns = split("\t", $line);
        ($header, $start, $end, $dir) = ($columns[0], ($columns[3] - 1), $columns[4], $columns[6]);
        $barrnap_28S{$header} = [$start, $end, $dir];
    }
}
        
close $in_barrnap;


# check blast and barrnap hash keys against each other and print values to a temp file
open(my $out_temp, ">$temp") or die "error creating $temp";

foreach my $header (keys %blast_coords) {
    if (exists $barrnap_18S{$header}) {
        print $out_temp "$header\t$barrnap_18S{$header}[2]\t$blast_coords{$header}[0]\t$blast_coords{$header}[1]\t$barrnap_18S{$header}[0]\t$barrnap_18S{$header}[1]\n";
    }
}

close $out_temp;


# open temp file and make hash that will be used to extract 18S sequences
open(my $in_temp, "<$temp") or die "error opening $temp";

while (my $line = <$in_temp>) {
    ($header, $dir, @columns) = split("\t", $line);
    my @ordered = (sort { $a <=> $b } @columns);
    my $min = $ordered[0];
    my $max = $ordered[-1];
    $extract{$header} = [$min, $max, $dir];
}

close $in_temp;


# open output file and extract 18S sequences using the substr function
open (my $output18S, ">$SSUout") or die "error opening $SSUout for writing";

foreach my $header (keys %extract) {
    if (exists $fasta{$header}) {
        print $output18S ">$header\n";
        if ($extract{$header}[2] eq "+") {
            my $seq = substr($fasta{$header}, 0, $extract{$header}[1]);
            print $output18S "$seq\n";
        } elsif ($extract{$header}[2] eq "-") { 
            my $seq = substr($fasta{$header}, 0, -$extract{$header}[0]);       
            print $output18S "$seq\n";
        }       
    }               
}


close $output18S;



# open output file and extract rrna sequences using the substr function
open (my $output28S, ">$LSUout") or die "error opening $LSUout for writing";

foreach my $header (keys %barrnap_28S) {
    if (exists $fasta{$header}) {
        print $output28S ">$header\n";
        if ($barrnap_28S{$header}[2] eq "-") {
            my $seq = substr($fasta{$header}, -$barrnap_28S{$header}[1]);
            print $output28S "$seq\n";
        } elsif ($barrnap_28S{$header}[2] eq "+") {
            my $seq = substr($fasta{$header}, $barrnap_28S{$header}[0]);
            print $output28S "$seq\n";
        }       
    }       
}

close $output28S;
