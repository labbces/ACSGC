#!/usr/bin/perl

#Requires MashMap/2.0 EMBOSS/6.0 gnuplot

use strict;
use warnings;
#use as:
#computeDotPlots.pl contigPairsFile.pl genomeFile.fasta
open TBL, $ARGV[0];
while(<TBL>){
 my ($contig1,undef,$contig2,undef,undef,undef,undef)=split(/\t/);
 my $fileContig1=$contig1.'.fasta';
 my $fileContig2=$contig2.'.fasta';
 my $BaseName=$contig1.'_vs_'.$contig2.'.mashmap';
 my $mashMapFile=$BaseName.'.out';
 my $mashMapFileLOG=$BaseName.'.out.log';
 my $dotPlotFilePS=$BaseName.'.dotplot.ps';
 my $dotPlotFileLOG=$BaseName.'.dotplot.log';
 if (!-f $fileContig1){
  system("seqret -auto $ARGV[1]:$contig1 $fileContig1")
 }
 if (!-f $fileContig2){
  system("seqret -auto $ARGV[1]:$contig2 $fileContig2")
 }
 system("mashmap -r $fileContig1 -q $fileContig2 --threads $ARGV[2] --output $mashMapFile -f none > $mashMapFileLOG 2>&1");
 system("generateDotPlot postscript large $mashMapFile > $dotPlotFileLOG 2>&1");
 if(-f 'out.ps'){
  rename 'out.ps', $dotPlotFilePS;
  rename 'out.rplot', $BaseName.'.rplot';
  rename 'out.fplot', $BaseName.'.fplot';
  rename 'out.gp', $BaseName.'.gp';
 }
}
