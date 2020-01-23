#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $debug=0;
my $clstrFile='';
my $gffFile='';
my $genomeFile='';
my $help='';
my $license='';
my $version='1.0.0';

############################################
##Get input from user
############################################

GetOptions ("gff|g=s"    => \$gffFile,
            "cluster|c=s"=> \$clstrFile,
            "genome|g=s" => \$genomeFile,
            "help|h|?"   => \$help, 
            "debug|d=i"  => \$debug,
            "license|l"  => \$license)
or die("Error in command line arguments\n");
############################################
##Check input from user
############################################
if($help){
 &usage;
 exit 1
}
if($license){
 &license;
 exit 1
}
if(!-s $gffFile){
 print STDERR "\n\tFATAL:  You must provide a GFF file with you genome annotation.\n\n";
 &usage;
 exit 0;
}
if(!-s $genomeFile){
 print STDERR "\n\tFATAL:  You must provide a GENOME file in fasta format. Make sure that the contigs/scaffold/chromosome identifiers match with those in your GFF file.\n\n";
 &usage;
 exit 0;
}
if(!$clstrFile){
 print STDERR "\n\tFATAL: You must provide the name of a file to store the results of running CD-hit on the predicted sets of proteins in your genome\n\n";
 &usage;
 exit 0;
}



#use like:
#./detectAllelicContigs.pl ../braker/augustus.hints.gff augustus.hints.nr90.aa.clstr
#

my %genes2contigs;
my %contigs2genes;
my %transcripts2clusters;
my %clusters2transcripts;
my %clusters2contigs;
my %contigs2clusters;

open GFF, $gffFile;
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

open CDHIT, $clstrFile;
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




############################################
#Usage
############################################
sub usage{
    print STDERR "$0 version $version, Copyright (C) 2015 Diego Mauricio Riano Pachon\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0 Detects putative pairs of allelic contigs based on shared gene content.

USAGE
    $0 --gff genome.gff3 --genome genome.fasta --cluster predicted_proteins.clstr

OPTIONS
    --help      -h    This help.
    --license   -l    License.

EOF
}
############################################
#License
############################################
sub license{
    print STDERR <<EOF;

Copyright (C) 2020 Diego Mauricio Riaño Pach<C3>
e-mail: diego.riano\@cena.usp.br

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
EOF
exit;
}

