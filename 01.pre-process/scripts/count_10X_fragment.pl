#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

### 10X fragment files
# chr1    3062785 3063037 TGTCTTGAGGCCCTTG-1      1
# chr1    3062786 3063059 TTCTTGCTCTTGTTGG-1      14

my %id;
my @exts = qw(.tsv .gz);
my($name, $path, $ext) = fileparse($ARGV[0], @exts);
if($ext eq ".gz"){
	open IN, "zcat $ARGV[0] |" or die $!;
}
else{
	open IN, "cat $ARGV[0] |" or die $!;
}

while(<IN>){
	next if $_ =~ m/^\#/;
	my @sp = split/\s+/, $_;
	my $cid = $sp[3];
	my $count = $sp[4];
	$id{$cid} = 0 if not exists $id{$cid};
	$id{$cid} = $id{$cid} + $count;
}
close IN;

open OUT, ">$ARGV[0]\_Count.xls" or die $!;
foreach my $i (keys %id){
	print OUT $i."\t$id{$i}\n";
}
close OUT;








