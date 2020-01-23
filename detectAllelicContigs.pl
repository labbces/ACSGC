#!/usr/bin/perl

use strict;
use warnings;

#use like:
#./detectAllelicContigs.pl ../braker/augustus.hints.gff augustus.hints.nr90.aa.clstr
#

my %genes2contigs;
my %contigs2genes;
my %transcripts2clusters;
my %clusters2transcripts;
my %clusters2contigs;
my %contigs2clusters;

open GFF, $ARGV[0];
while(<GFF>){
 chomp;
 next if /^#/;
 my @fields=split("\t");
 if ($fields[2] eq 'transcript'){
  $genes2contigs{$fields[8]}=$fields[0];
  $contigs2genes{$fields[0]}{$fields[8]}=1;
 }
}
close GFF;

open CDHIT, $ARGV[1];
my $clusterId='';
while(<CDHIT>){
 chomp;
 if (/^>Cluster (\d+)$/){
  $clusterId='cluster_'.$1;
 }
 else{
  my @f1=split(/,/);
  my @f2=split(/ /,$f1[1]);
  my $seqId=$f2[1];
  $seqId=~s/^>//;
  $seqId=~s/\.+$//;
  $transcripts2clusters{$seqId}=$clusterId;
  $clusters2transcripts{$clusterId}{$seqId}=1;
 }
}
close CDHIT;

foreach my $cluster(keys %clusters2transcripts){
 foreach my $seqID(keys(%{$clusters2transcripts{$cluster}})){
  $clusters2contigs{$cluster}{$genes2contigs{$seqID}}=1;
  $contigs2clusters{$genes2contigs{$seqID}}{$cluster}=1;
 }  
}

my %selectedContigs;
foreach my $cluster(keys %clusters2contigs){
 if(keys(%{$clusters2contigs{$cluster}}) > 1){
  foreach my $contig(keys(%{$clusters2contigs{$cluster}})){
   $selectedContigs{$contig}=1;
  }
 }
}

print STDERR "There are ". scalar(keys(%selectedContigs))." contigs with shared gene/protein clusters\n";
foreach my $contig1(keys %selectedContigs){
# my @clustersContigs1=keys %{$contigs2clusters{$contig1}};
 my @clustersContigs1=keys %{$contigs2clusters{$contig1}};
 my %clustersContigs1 = map { $clustersContigs1[$_] => 1 } 0 .. $#clustersContigs1;
 foreach my $contig2(keys %selectedContigs){
  my $countShared=0;
  next if $contig1 eq $contig2;
  my @clustersContigs2=keys %{$contigs2clusters{$contig2}};
  my %clustersContigs2 = map { $clustersContigs2[$_] => 1 } 0 .. $#clustersContigs2;
  foreach my $cc1(keys %clustersContigs1){
   if($clustersContigs2{$cc1}){
    $countShared++;
   }
  }
  my $fractionSharedContig1=($countShared*100)/scalar(keys(%{$contigs2clusters{$contig1}}));
  my $fractionSharedContig2=($countShared*100)/scalar(keys(%{$contigs2clusters{$contig2}}));
  #Selecting pairs of contigs if at least one of the contigs shares 50% of their genes with the other contig
  if(($fractionSharedContig1 >=50 || $fractionSharedContig2 >=50) && (@clustersContigs1 > 10 && @clustersContigs2 > 10)){
   print "$contig1\t".scalar(keys(%{$contigs2clusters{$contig1}}))."\t$contig2\t".scalar(keys(%{$contigs2clusters{$contig2}}))."\t$countShared\t$fractionSharedContig1\t$fractionSharedContig2\n";
  }
 }
}
