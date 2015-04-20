import argparse
import pysam
import multiprocessing as mp
import signal

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-v', '--variants', required=True)
parser.add_argument('-r', '--reads', required=True)
parser.add_argument('-p', '--processes', type=int, default=4)
parser.add_argument('--qscore-threshold', type=int, default=15)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# determine variant positions

def main():
  # load all the regions
  bamfile = pysam.Samfile(args.bam, "rb")
  contigs = [contig for contig in bamfile.references]

  # launch subprocesses
  pool = mp.Pool(processes=args.processes)
  try:
    results = pool.map(parse_region, contigs)
  except KeyboardInterrupt:
    pool.terminate()
    exit("Pool terminated by user")

  with open(args.reads, 'w') as out:
    for contig, reads in zip(contigs, results):
      for read in reads:
        out.write('%s\t%s' % (contig, read))
        for pos, seq in sorted(reads[read].iteritems()):
          out.write('\t%d:%s' % (pos, seq))
        out.write('\n')

def parse_region(contig):
  signal.signal(signal.SIGINT, signal.SIG_IGN)

  # load all variants into shared memory
  n = 0
  variants = dict()
  with open(args.variants) as f:
    for line in f:
      fields = line.strip().split()
      ctg, pos, alleles = fields[0], int(fields[1]), {fi.split(':')[0] for fi in fields[2:]}

      if ctg != contig:
        continue
      else:
        variants[pos] = alleles

  bamfile = pysam.Samfile(args.bam, "rb")
  my_reads = dict()
  for pos in variants:
    for pcol in bamfile.pileup(contig, start=pos, end=pos+1, truncate=True):
      alleles = variants[pcol.pos]

      for pread in pcol.pileups:
        if pread.alignment.query_qualities[pread.query_position] < args.qscore_threshold:
          continue
        if pread.is_del:
          # deletion
          seq = "-"
        elif pread.indel > 0:
          # we have an insertion
          seq = 'I' + pread.alignment.seq[pread.query_position:pread.query_position+1+pread.indel]
        else:
          # normal match
          seq = pread.alignment.seq[pread.query_position]

        if seq in alleles:
          read = pread.alignment.qname
          if read not in my_reads:
            my_reads[read] = dict()
          my_reads[read][pcol.pos] = seq

  return my_reads

if __name__ == '__main__':
  main()
