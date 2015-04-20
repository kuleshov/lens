import argparse
import pysam
import signal
import multiprocessing as mp

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-r', '--reads', required=True)
parser.add_argument('-k', '--clusters', required=True)
parser.add_argument('-p', '--processes', type=int, default=4)
parser.add_argument('--cov-cutoff', default=2, type=int,
                    help='minimum coverage to call a haplotype')
parser.add_argument('--cov-cutoff-percentage', default=0.75, type=float,
                    help='percentage of haplotype that needs to be covered with cov-cutoff')
parser.add_argument('--similarity-cutoff', default=1.0, type=float,
                    help='minimum % similarity to combine two haplotypes')
parser.add_argument('--overlap-cutoff', default=2, type=int,
                    help='minimum number of overlaping positions necessary to combine two haplotypes')

args = parser.parse_args()

# ----------------------------------------------------------------------------
# helper classes/functions

class Profile:
  def __init__(self, profile):
    # _profile is a dict of the form {pos:allele}
    assert profile
    self._profile = profile

  def common_pos(self, other):
    return [allele for allele in self._profile if allele in other._profile]

  def similar(self, profile):
    n_equal = len([x for x in self if x in profile])
    n_common = len(self.common_pos(profile))
    if n_common >= args.overlap_cutoff and (n_equal >= args.similarity_cutoff*n_common):
      return True
    else:
      return False

  def __len__(self):
    return len(self._profile)

  def __iter__(self):
    return self._profile.iteritems()

  def __contains__(self, x):
    return (self._profile.get(x[0], None) == x[1])

  def __eq__(self, other):
    fz0 = frozenset(self._profile.items())
    fz1 = frozenset(other._profile.items())
    return fz0 == fz1

  def __hash__(self):
    return frozenset(self._profile.items()).__hash__()

  def update(self, profile):
    new_profile = dict(self._profile.items() + profile._profile.items())
    return Profile(new_profile)

  def start(self):
    return sorted(self._profile.keys())[0]

  def end(self):
    return sorted(self._profile.keys(), reverse=True)[0]

  def __str__(self):
    sorted_items = sorted((x for x in self._profile.iteritems()), key=lambda x: x[0])
    return ''.join(x[1] for x in sorted_items)

  def string_over(self, positions):
    return ''.join([self._profile.get(pos, '?') for pos in sorted(positions)])

  def delete(self, pos):
    del self._profile[pos]

def write_clusters(cluster_file, clusters, contig, reads, n):
  sorted_items = sorted(clusters.iteritems(), key=lambda x: (x[0].start(), x[0].end()))
  region_positions = {p for profile, R in sorted_items for (p,x) in profile}
  for profile, R in sorted_items:
    if len(R) < 3:
      continue

    block_positions = sorted(profile._profile.keys())
    coverages = {pos:0 for pos in block_positions}
    for read in R:
      for pos in block_positions:
        if pos in reads[read]: coverages[pos] += 1
    coverage = float(sum(coverages.values())) / len(coverages)

    for read in R:
      read_seq = [reads[read].get(pos, '?') for pos in sorted(profile._profile.keys())]

    if (len([c for c in coverages.values() if c >= args.cov_cutoff]) < args.cov_cutoff_percentage*len(block_positions)):
      continue

    for pos in block_positions:
      if coverages[pos] < 2:
        profile.delete(pos)

    cluster_name = contig + '-' + str(n)
    cluster_file.write('%s\t%d\t%d\t%f\t' % (contig, profile.start(), profile.end(), coverage))
    pos_string = ','.join([('%d:%s' % (p,a)) for (p,a) in sorted(profile)])
    cluster_file.write(pos_string)
    for read in R:
      cluster_file.write('\t%s' % read)
      read_seq = [reads[read].get(pos, '?') for pos in sorted(profile._profile.keys())]

    #print
    cluster_file.write('\n')

# ----------------------------------------------------------------------------
# determine groups of reads

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

  # load reads
  reads = dict()
  with open(args.reads) as f:
    for line in f:
      fields = line.strip().split()
      ctg, read = fields[0], fields[1]
      if ctg not in reads:
        reads[ctg] = dict()
      alleles = {int(fi.split(':')[0]):fi.split(':')[1] for fi in fields[2:]}
      reads[ctg][read] = alleles

  with open(args.clusters, 'w') as out:
    for contig, cluster_groups in zip(contigs, results):
      for n, clusters in enumerate(cluster_groups):
        write_clusters(out, clusters, contig, reads[contig], n)

def parse_region(contig):
  signal.signal(signal.SIGINT, signal.SIG_IGN)
  reads = dict()
  with open(args.reads) as f:
    for line in f:
      fields = line.strip().split()
      ctg, read = fields[0], fields[1]
      if ctg != contig:
        continue

      alleles = {int(fi.split(':')[0]):fi.split(':')[1] for fi in fields[2:]}
      reads[read] = alleles
  bamfile = pysam.Samfile(args.bam, 'rb')

  curr_reads = set()
  cluster_groups = list()
  clusters = dict()

  poslist = {p for r in reads for p in reads[r]}
  n = 0
  # move from left to right and add reads to current group
  for pos in sorted(poslist):
    # fetch new reads at this position:
    reads_at_pos = {r.qname for r in bamfile.fetch(contig, pos, pos+1)}
    old_reads = reads_at_pos & curr_reads
    new_reads = reads_at_pos - curr_reads

    if curr_reads and not old_reads:
      cluster_groups.append(clusters)
      n += 1
      clusters = dict()
      curr_reads = set()

    # assign new reads to clusters, or form new clusters
    for r in new_reads:
      if r in reads:
        if not reads[r]: continue
        read_profile = Profile(reads[r])
      else:
        continue
      for cluster_profile in clusters:
        if read_profile.similar(cluster_profile):
          #print read_profile, ' similar to ', cluster_profile
          new_profile = cluster_profile.update(read_profile)
          clusters[new_profile] = clusters.pop(cluster_profile)
          clusters[new_profile].add(r)
          break
      else:
        clusters[read_profile] = set([r])

    curr_reads.update(new_reads)

  if clusters:
    cluster_groups.append(clusters)

  return cluster_groups

if __name__ == '__main__':
  main()
