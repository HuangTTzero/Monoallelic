#!/usr/bin/perl
use warnings;
use strict;

### This script filters the latest VCF file for various SNPs and writes high confidence SNPs into a folder called 'SNPs_Sanger' (needs to exist already)
### Here we are interested in the strain 129S1_SvImJ for a collaboration with Andrew Dimond

### last modified 11 April 2014

#ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz
my $infile = '/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz'; 
warn "Processing SNPs from $infile\n\n";

if ($infile =~ /gz$/){
  open (IN,"zcat $infile |") or die "Failed to open file '$infile': $!\n";
}
else{
  open IN, $infile or die $!;
}

my $count = 0;

my $other = 0;
my $too_many = 0;

my %fhs;
my $hcg_count = 0;
my $low_confidence = 0;
my $same = 0;
my $homozygous = 0;
my %all_SNPs;
my @origin_SNPs=();

for my $chr (1..19,'X','Y') {
  my $filename = 'SNPs_Sanger/chr'.$chr.'.txt';
  open (my $fh,'>',$filename) or die "Couldn't open filehandle $!\n";
  $fhs{$chr} = $fh;
  print {$fhs{$chr}} ">$chr\n";
}

#for (1..11){
#  $_ = <IN>; #headers and info lines
#}

while (<IN>){
  chomp;
  next if ($_ =~ /^\#\#/);
  $count++;
  if ($count%1000000 ==0){
    warn "processed $count lines\n";
  }
  #  last if ($count == 1000000);
  #  warn "$_\n"; sleep(1);
  my ($chr,$pos,$ref,$alt,$strain129) = (split /\t/)[0,1,3,4,9];
  # warn "$strain129\n";
  # sleep(5);

  my ($gt,$gq,$dp,$pl,$sp,$fi) = (split/:/,$strain129)[0,1,2,5,10,13];
  unless ($fi){
    # warn "genotype: $gt\nfilter:   $fi\n\n";
  }
  # $gt is the Genotype:


  # '.'   = no genotype call was made
  # '0/0' = genotype is the same as the reference geneome
  # '1/1' = homozygous alternative allele; can also be '2/2',
  # '3/3', etc. if more than one alternative allele is present.
  # '0/1' = heterozygous genotype; can also be '1/2', '0/2', etc.

  # $fi is filter, 1 for high confidence SNP, or 0 for low confidence

  ### I guess we are only looking for 1/1 calls, and filter for high confidence as well

  if ($ref =~ /[^ATCG]/){
    warn "ref was: $ref; skipping\n";
    next;
  }

  # skipping if the SNP is not well defined in C3H
  if ($alt =~ /[^ATCG]/){
    # warn "SNP was: $alt; skipping; $c3h\n";
    ++$too_many;
    next;
  }

  if ($gt eq '0/0'){
    ++$same;
    # warn "same as reference\n";
    next;
  }
  elsif ($gt eq '1/1'){
    ++$homozygous;
    # warn "homozygous alternative allele\n";
  }
  else{
    ++$other;
    next;
  }

  # Looking at the Filtering tag

  if ($fi == 1){
    ++$hcg_count;

    my $location = join (':',$chr,$pos);
    # warn "$location\n";
    my $SNP = join ("\t",$count,$chr,$pos,'1',join ("\/",$ref,$alt));
    # warn "$SNP\n";sleep(1);

    if (exists $all_SNPs{$location} ){
      warn "meep\n";
    }
    else{
      $all_SNPs{$location} = $SNP;
      push @origin_SNPs,$_;
    }


  }
  else{
    ++$low_confidence;
    next;
  }

  # Output example
  # Variation ID    Chromosome name Position on Chromosome (bp)     Strand  Allele
  # rs2020560       10      98212004        1       A/T
  print {$fhs{$chr}} join ("\t",$count,$chr,$pos,'1',join ("\/",$ref,$alt),$strain129),"\n";
}

warn "\n\n=========================================================\n\nPositions read in total:\t$count\n\n";
warn "SNP position summary for strain PWK\n\n";

warn "$homozygous\tSNP were homozygous. Of these:\n";
warn "$hcg_count\tSNP were homozygous and passed high confidence filters, and were included into PWK genome\n";

warn "\nNot included into PWK genome:\n";
warn "$same\thad the same sequence as the reference\n";
warn "$too_many\thad no clearly defined alternative base\n";
warn "$other\tCalls were neither 0/0 (same as reference) or 1/1 (homozygous SNP)\n";
warn "$low_confidence\twere homozygous but the filtering call was low confidence\n";


warn "Now printing a single list of all SNPs\n";
open (OUT,'>','all_PWK_SNPs.txt') or die $!;

foreach my $location (keys %all_SNPs){
  print OUT "$all_SNPs{$location}\n";
}
open (OUT,'>','origin_PWK_SNPs.txt') or die $!;

foreach my $rawline(@origin_SNPs){
  print OUT "$rawline\n";
}
warn "Done\n\n";
