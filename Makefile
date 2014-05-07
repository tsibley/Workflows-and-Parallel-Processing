SHELL:=/bin/bash
export SHELLOPTS:=errexit:pipefail
export PATH:=.:$(PATH)
.DELETE_ON_ERROR:

%_aa.fa: %_na.fa
	transeq -sequence $< -outseq $@ \
			-frame 1 -clean

%.nxs: %.fa
	muscle -quiet < $< | fasta2nexus > $@

# Keep intermediate alignments, for speed
.PRECIOUS: %.nxs

%_aa_freq.tsv: %_aa.nxs
	perl CountAAFreq.pl $< $@ 0.25 0.5
