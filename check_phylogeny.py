#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-k', '--clusters', required=True)
parser.add_argument('-p', '--perfect', required=True)
parser.add_argument('-s', '--stats', required=True)
parser.add_argument('-o', '--pos-stats', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# extract clusters

clusters = dict()
contigs = dict()
with open(args.clusters) as f:
  for line in f:
    fields = line.split()
    id = fields[3]
    positions = {int(fi.split(':')[0]):fi.split(':')[1] for fi in fields[6].split(',') if fi}
    if id not in clusters: clusters[id] = list()
    clusters[id].append(positions)
    contigs[id] = fields[0]

# ----------------------------------------------------------------------------
# functions

def write_cluster(name, cluster, all_pos, out):
  start = min(all_pos)
  end = max(all_pos)
  for D in cluster:
    string = ''.join([D.get(pos, '-') for pos in all_pos])
    out.write('%s\t%d\t%d\t%s\n' % (name, start, end, string))

# ----------------------------------------------------------------------------
# count number of perfect phylogenies

discarded = list()
perfect = list()
imperfect = list()

perfect_file = open(args.perfect, 'w')
pos_stats = open(args.pos_stats, 'w')

for name, cluster in clusters.iteritems():
  # extract all positions and alleles
  all_pos = dict()
  for D in cluster:
    for pos, a in D.iteritems():
      if pos not in all_pos: all_pos[pos] = set()
      all_pos[pos].add(a)

  # common_pos = all_pos.keys()
  # for D in cluster:
  #   common_pos = 

  counts = {pos:len(L) for (pos, L) in all_pos.iteritems()}
  bad_pos = [pos for pos in all_pos if counts[pos] > 2]
  for pos in bad_pos:
    pos_stats.write('%s\t%d\n' % (contigs[name], pos))


  if max(counts.values()) > 2 and float(len(bad_pos)) / len(counts) > 0.05:
    # discard clusters with more than 2 alleles at a position
    # print name, len(), len(counts)
    discarded.append((len(all_pos), len(cluster)))
    continue
  else:
    # translate alleles to 0,1
    tr = {pos:{a:i for (i,a) in enumerate(all_pos[pos])} for pos in all_pos if pos not in bad_pos}
    good_pos = [pos for pos in all_pos if pos not in bad_pos]

  # build matrix
  M = list()
  for i, D in enumerate(cluster):
    M.append(list())
    for pos in sorted(good_pos):
      if pos in D:
        M[i].append(tr[pos][D[pos]])
      else:
        M[i].append('-')

  # does it satisfy phylogeny?
  sat = True
  problem_pairs = set()
  for coli in xrange(len(M[0])):
    for colj in xrange(len(M[0])):
      if coli < colj: continue
      tuples = set([(M[k][coli], M[k][colj]) for k in xrange(len(cluster))])
      if (0,1) in tuples and (1,0) in tuples and (0,0) in tuples and (1,1) in tuples:
        sat = False
        problem_pairs.add((coli, colj))

  if sat:
    perfect.append((len(all_pos), len(cluster)))
    write_cluster(name, cluster, all_pos.keys(), perfect_file)
  else:
    imperfect.append((len(all_pos), len(cluster), len(problem_pairs)))
    sorted_pos = sorted(tr.keys())
    prob_cols = set([c for p in problem_pairs for c in p])
    for c in prob_cols: pos_stats.write('%s\t%d\n' % (contigs[name], pos))

with open(args.stats,'w') as stats:
  stats.write('# counts of % of wrong column pairs:\n')
  for x in perfect:
    stats.write('0 ')
  for x in imperfect:
    stats.write('%f ' % ( float(x[2]) / x[0] ) )
  stats.write('\n')

  stats.write('# correct and incorrect (>2% of pairs wrong) pts and the length and depth:\n')
  for x in perfect:
    stats.write('0.0,%d,%d ' % x)
  for x in imperfect:
    if ( float(x[2]) / x[0] <= 0.02):
      stats.write('%f,%d,%d ' % (float(x[2]) / x[0], x[0], x[1]))
  stats.write('\n')
  for x in imperfect:
    if ( float(x[2]) / (x[0]) > 0.02):
      stats.write('%f,%d,%d ' % (float(x[2]) / x[0], x[0],x[1]))
  stats.write('\n')
  for x in discarded:
    stats.write('-1.0,%d,%d ' % (x))
  stats.write('\n')

print len(perfect), len(imperfect), len(discarded)
print sum([x[0] for x in perfect])/len(perfect), sum([x[1] for x in perfect])/len(perfect), max([x[1] for x in perfect])
print sum([x[0] for x in imperfect])/len(imperfect), sum([x[1] for x in imperfect])/len(imperfect), max([x[1] for x in imperfect])
print sum([x[0] for x in discarded])/len(discarded), sum([x[1] for x in discarded])/len(discarded), max([x[1] for x in discarded])
print [x for x in imperfect]

