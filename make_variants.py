import argparse
import pysam
import multiprocessing as mp
import signal

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-o', '--pos', required=True)
parser.add_argument('-p', '--processes', type=int, default=4)
parser.add_argument('--indels', action='store_true',
                    help='Call indels')
parser.add_argument('--coverage-threshold', type=int, default=3,
                    help='Minimum number of reads that need to support a variant')
parser.add_argument('--frequency-threshold', type=float, default=0.1,
                    help='Minimum allele frequency of a variant')
parser.add_argument('--qscore-threshold', type=int, default=15,
                    help='Minimum qscore at a position in a read for it to be considered')

args = parser.parse_args()

# ----------------------------------------------------------------------------
# determine variant positions

def main():
  bamfile = pysam.Samfile(args.bam, "rb")

  regions = list()
  for contig, length in zip(bamfile.references, bamfile.lengths):
    regions.extend([(contig, N+1, N+10000) for N in range(0,length,10000)])

  pool = mp.Pool(processes=args.processes)

  try:
    results = pool.map(parse_region, regions)
  except KeyboardInterrupt:
    pool.terminate()
    exit("Pool terminated by user")

  posfile = open(args.pos, "w")
  for (contig, start, end), results in zip(regions, results):
    for pos, significant_variants in sorted(results.iteritems()):
      posfile.write('%s\t%d' % (contig, pos))
      for b, n in sorted(significant_variants.iteritems(), key=lambda x: x[1], reverse=True):
        posfile.write('\t%s:%d' % (b, n))
      posfile.write('\n')

  posfile.close()


def parse_region(arg):
  signal.signal(signal.SIGINT, signal.SIG_IGN)

  contig, start, end = arg
  bamfile = pysam.Samfile(args.bam, "rb")
  variant_calls = dict()
  pileup = bamfile.pileup(reference=contig, start=start, end=end, truncate=True)
  for pcol in pileup:
    alleles = dict()
    for pread in pcol.pileups:
      if pread.alignment.query_qualities[pread.query_position] < args.qscore_threshold: continue
      if not args.indels:
        if pread.indel != 0: continue
        if pread.is_del: continue
      if pread.is_del:
        # deletion
        seq = "-"
      elif pread.indel > 0:
        # we have an insertion
        seq = 'I' + pread.alignment.seq[pread.query_position:pread.query_position+1+pread.indel]
      else:
        # normal match
        seq = pread.alignment.query_sequence[pread.query_position]

      if seq not in ('N',):
        if seq not in alleles:
          alleles[seq] = 0
        alleles[seq] += 1

    total_alleles = sum(alleles.values())
    threshold = max(args.coverage_threshold, args.frequency_threshold*total_alleles)
    significant_variants = {k:v for (k,v) in alleles.iteritems() if v >= threshold}
    if len(significant_variants) >= 2:
      variant_calls[pcol.pos] = significant_variants

  return variant_calls

if __name__ == '__main__':
  main()



