# USES THE LENS ALGORITHM TO PRODUCE BACTERIAL HAPLOTYPES WITHIN CONTIGS
# Input files: sample.bam
# Output files: haplotypes.txt

#------------------------------------------------------------------------------
# High-level commands

run:
	make haplotypes.txt

clean:
	-rm variants.* haplotypes.txt

#------------------------------------------------------------------------------
# Alignment process

BAM:=./sample.bam
PROCESSORS:=8
LENS:=..

variants.pos: sample.bam sample.bam.bai
	python $(LENS)/make_variants.py -b $< -p $(PROCESSORS) -o $@.tmp
	mv $@.tmp $@

variants.reads: sample.bam sample.bam.bai variants.pos
	python $(LENS)/make_reads.py -b $< -p $(PROCESSORS) -v variants.pos -r $@.tmp
	mv $@.tmp $@

haplotypes.txt: sample.bam variants.pos variants.reads
	python $(LENS)/detect_subspecies.py -b $< -p $(PROCESSORS) -r variants.reads -k $@.tmp
	mv $@.tmp $@
