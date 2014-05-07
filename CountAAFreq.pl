#!/usr/bin/perl -w

#####################################################################################################
# Program: CountAAFreq.pl
# Purpose: count the frequency of amino acid at each position of sequence in nexus alignment file
# Input: MacClade nexus alignment file
# Output: tab delimited list of amino acid frequency at each position
# Author: Wenjie Deng
# Date: 2006-05-03
# Modified: 2006-05-15
# Changes: handle the alignment that there is no amino acids in entire columns, and 
# shift the column and output the correct position
# Modified: 2011-07-13
# Changes: handle sequences in alignment with different start and/or end positions
# Usage: perl CountAAFreq.pl inputFile outFile gapCutoff frequencyCutoff
#####################################################################################################

use strict;

my $usage = "Usage: perl countAaFreq.pl inputFile outFile gapCutoff(0-1) frequencyCutoff(0-1)\n";
my $inFile = shift or die $usage;	# amino acid sequence nexus file
my $unixInFile = $inFile."_unix";
my $outFile = shift or die $usage;	# output tab delimited list of amino acid frequency at each position
my $gapCutoff = shift or die $usage;
my $freqCutoff = shift or die $usage;
my $seqArr = ();
my @aaArr = ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-");
my $ntax = my $nchar = my $element = my $count = my $flag = 0;
my ($naStart, $naEnd);
system("tr '\r' '\n' < $inFile > $unixInFile");
open (IN, $unixInFile) or die "Couldn't open $unixInFile: $!\n";
open (OUT, ">$outFile") or die "Couldn't open $outFile: $!\n";

while (my $line = <IN>) {
	chomp $line;
	if ($line) {
		if ($line =~ /NTAX=(\d+)/i) {
			$ntax = $1;
			print "ntax: $ntax\n";
		}
		if ($line =~ /NCHAR=(\d+)/i) {			
			$nchar = $1;
			print "nchar: $nchar\n";
		}
		if ($line =~ /^\s*MATRIX\s*$/i) {
			$flag = 1;
		}elsif ($line eq ';') {
			$flag = 0;
		}elsif ($flag) {
			unless ($line =~ /^\[\s+/) {
				if ($line =~ /^['"]?(.*)['"]\s+(.*)(?:\s+\[\d+\])?$/) {
					my $aaSeq = uc $2;
					$aaSeq =~ s/\s//g;
					for (my $i = 0; $i < $nchar; $i++) {
						$seqArr->[$i]->[$element] = substr($aaSeq, $i, 1);
					}
					for (my $i = 0; $i < $nchar; $i++) {
						if ($seqArr->[$i]->[$element] =~ /[A-Za-z]/) {
							$naStart->{$element} = $i;
							last;
						}
					}
					for (my $i = $nchar - 1; $i >= 0; $i--) {
						if ($seqArr->[$i]->[$element] =~ /[A-Za-z]/) {
							$naEnd->{$element} = $i;
							last;
						}
					}	
					$element++;	# sequence index
				}
			}
		}	
	}	
}
close IN;

my $freqCutoffPercent = $freqCutoff * 100;
print OUT "Position\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\t-\t";
print OUT "Total\t\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\t\tMostFreqAa\tFrequency\t2nd AA\tFrequency\t$freqCutoffPercent","% AA\t$freqCutoffPercent","% Freq.\n";
my $position = 1;
for (my $i = 0; $i < $nchar; $i++) {
	my %aaHash = my %freqHash = my @mostFreqAAs = ();
	my $totalCount = my $flag = 0;	
	foreach my $aa (@aaArr) {
		$count = 0;
		$flag++;
		for (my $j = 0; $j < $ntax; $j++) {
			if ($i >= $naStart->{$j} && $i <= $naEnd->{$j}) {
				if ($flag == 1) {
					$totalCount++;
				}				
				if ($seqArr->[$i]->[$j] eq $aa) {
					$count++;
				}
			}
		}
		$aaHash{$aa} = $count;
	}
	my $totalAA = $totalCount - $aaHash{"-"};
	unless ($totalAA == 0) {
		print OUT $position,"\t";
		foreach my $aa (@aaArr) {
			print OUT $aaHash{$aa},"\t";
		}
		print OUT "$totalAA\t";
		my %aaFreqHash;
		my $topFreq = 0;
		my $topFreqAa = "";
		foreach my $aa (@aaArr) {
			unless ($aa eq "-") {
				my $frequency = int (10000 * $aaHash{$aa}/$totalAA + 0.5) / 10000;
				print OUT "\t$frequency";
				$freqHash{$aa} = $frequency;
			}
		}
		my $gapFreq = int (10000 * $aaHash{"-"}/$totalCount + 0.5) / 10000;
		if ($gapFreq > $gapCutoff) {
			$topFreqAa = "-";
			$topFreq = 0;
			print OUT "\t\t$topFreqAa\t$topFreq\n";
		}else {
			@mostFreqAAs = sort {$freqHash{$b} <=> $freqHash{$a}} keys %freqHash;
			if ($freqHash{$mostFreqAAs[0]} == 1) {
				print OUT "\t\t$mostFreqAAs[0]\t$freqHash{$mostFreqAAs[0]}\t\t\t$mostFreqAAs[0]\t$freqHash{$mostFreqAAs[0]}\n";
			}elsif ($freqHash{$mostFreqAAs[0]} >= $freqCutoff) {
				print OUT "\t\t$mostFreqAAs[0]\t$freqHash{$mostFreqAAs[0]}\t$mostFreqAAs[1]\t$freqHash{$mostFreqAAs[1]}\t$mostFreqAAs[0]\t$freqHash{$mostFreqAAs[0]}\n";
			}else {
				my $bothAA = $mostFreqAAs[0].", ".$mostFreqAAs[1];
				my $bothFreq = $freqHash{$mostFreqAAs[0]} + $freqHash{$mostFreqAAs[1]};
				if ($bothFreq < $freqCutoff) {
					$bothAA = $bothFreq = "";
				}
				print OUT "\t\t$mostFreqAAs[0]\t$freqHash{$mostFreqAAs[0]}\t$mostFreqAAs[1]\t$freqHash{$mostFreqAAs[1]}\t$bothAA\t$bothFreq\n";
			}
		}		
		$position++;
	}
}
close OUT;
unlink ($unixInFile);

print "Inplemented $ntax sequences, each sequence contained $nchar characters. All done!\n";
