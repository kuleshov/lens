Lens
====

Algorithm for uncovering genomic variants and haplotypes in metagenomes from long reads. Lens was originally introduced in

```
High-resolution structure of the human microbiome revealed with synthetic long reads.
Volodymyr Kuleshov, Chao Jiang, Wenyu Zhou, Fereshteh Jahanbani, Serafim Batzoglou, Michael Snyder. 
Nature Biotechnology, 2015.
```

Rich bacterial communities, such as the human microbiome,
are composed of thousands of stains;
these strains may greatly differ from each other only
in their virulence or in their resistance to antibiotics.

Modern sequencing technologies have difficulty finding
differences between closely related strains at the SNP 
resolution. This problem however becomes much easier with longer
reads such as ones obtained from recent synthetic long read 
library preparation protocols.

Lens is a software package that detects fine-grained variation
in metagenomic samples based on long read data and
assembles these variants into long haplotypes.

At the core of Lens is a greedy phasing algorithm that
generalizes single-individual haplotyping algorithms to
the setting where there the number of underlying haplotypes
may be greater than two, and moreover may be unknown.

## Installation

To install Lens, simply clone the git repo:

```
git clone [...];
cd lens;
```

### Requirements

Lens requires `python-2.7` and the package `pysam` version `>=0.8.2`. 
The easiest way to install it using `pip install pysam`.

### Testing the installation

To test whether Lens works, we have provided a small testing package in 
`lens/test`:

```
cd lens/test;
make run;
```

This executes Lens on a small `.bam` containing real synthetic long reads.
The reults is a file called `haplotypes.txt`, which should contain 6 haplotypes
recovered from the long read data.

## Running Lens

Lens takes as input a `.bam` file with long reads aligned to a reference genome
or to an assembled genomic contig. We recommend using `bwa mem` for generating this
file.

###### Step 0: BAM filtering

We suggest running Lens only on reads with high mapping scores (`>=30`).
We also suggest using reads that align across their entire length (that aren't clipped).
We filter for such reads using:

```
python filter_by_cigar.py -i input.bam -o filtered.bam;
samtools view -b -h -q 30 filtered.bam > filtered.q30.bam
```

###### Step 1: Variant calling

Our first step is to call variants in the reference:

```
python make_variants.py \
  -b filtered.q30.bam \
  -o variants.pos \
  --coverage-threshold 3 \
  --frequency-threshold 0.1 \
  --qscore-threshold 15 \
  --indels
```

This will find all variants (including indels) that have at least
3 reads supporting them and have an allele frequency of at least 0.1
(these are the default parameters). Also, only reads that have a
qscore of 15 or more at a given position will be used for 
calculating support.

###### Step 2: Variant calling in reads

Next, we compile the reads that support the variants we have just  identified.
The result of this step is a file that lists for each read all the variants
that the read covers.

```
python make_reads.py \
  -b filtered.q30.bam \
  -v variants.pos \
  -r variants.reads \
  --qscore-threshold 15
```

The format of `variants.reads` is:
```
<ctg id>   <read id>   <pos1:allele1>   <pos2:allele2>  ...
```

###### Step 3: Assembling reads into haplotypes

Finally, we assemble the reads at their overlapping variants into bacterial haplotypes.

```
python detect_subspecies.py \
  -b filtered.q30.bam \
  -r variants.reads \
  -k haplotypes.txt \
  --cov-cutoff 2 \
  --cov-cutoff-percentage 0.75 \
  --similarity-cutoff 1.0 \
  --overlap-cutoff 2
```

This will report only haplotypes that have a coverage of at least two over at least
75% of their length. Two reads will be connected into a haplotype if they overlap 
at at least two variants and are identical at both of these positions.

The format of `haplotypes.txt` is:
```
<ctg id>   <start pos>   <end pos>   <coverage>   <variants>
```

## Feedback

Please send feedback to [Volodymyr Kuleshov](http://www.stanford.edu/~kuleshov)
