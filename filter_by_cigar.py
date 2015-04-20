import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input-bam', required=True)
parser.add_argument('-o', '--output-bam', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# code

infile = pysam.Samfile(args.input_bam, "rb")
outfile = pysam.Samfile(args.output_bam, "wb", template=infile)

n0 = 0
n1 = 0
no_cigar = 0
for read in infile:
  if not read.cigar:
    no_cigar += 1
    continue

  clipped_bases = sum([l for (o,l) in read.cigar if o in (4, 5)])
  if clipped_bases < 500:
    n1 += 1
    outfile.write(read)

  n0 += 1
  # if n0 % 10000 == 0: print n0

if no_cigar:
  print "WARNING: %d reads did not have a cigar string" % no_cigar

print "%d / %d reads passed the CIGAR filter" % (n1, n0)

infile.close()
outfile.close()
