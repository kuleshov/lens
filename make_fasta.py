#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--clusters', required=True)
parser.add_argument('-f', '--fasta', required=True)
parser.add_argument('-o', '--out-fasta', required=True)
parser.add_argument('-t', '--cluster-variants', required=True)
parser.add_argument('--variant-subset')
parser.add_argument('--indels', action='store_true')
# parser.add_argument('-c', '--contig', required=True)
# parser.add_argument('-s', '--start', required=True, type=int)
# parser.add_argument('-e', '--end', required=True, type=int)
# parser.add_argument('--snps', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

fasta = pysam.FastaFile(args.fasta)
out = open(args.out_fasta, 'w')
out_variants = open(args.cluster_variants, 'w')

variant_subset = set()
if args.variant_subset:
  with open (args.variant_subset) as f:
    for line in f:
      fields = line.strip().split()
      variant_subset.add((fields[0], int(fields[1])))

with open(args.clusters) as f:
  n=0
  m=0
  for line in f:
    fields = line.split()
    if len(fields) <= 5: continue
    contig, start0, end0 = fields[0], int(fields[1]), int(fields[2])
    # start, end = max(start0-75,0), end0+75
    start, end = start0, end0
    fields[4] = fields[4].strip(',')
    snps = dict([f.split(':') for f in fields[4].split(',')])
    seq = list(fasta.fetch(contig, start, end+1))
    tr_snps = dict()
    delta = 0
    for pos, snp in sorted(snps.iteritems()):
      delta_pos = 0
      if (snp == '-' or snp.startswith('I')): 
        is_indel = True
      else:
        is_indel = False

      if is_indel and not args.indels:
        continue

      pos_in_seq = int(pos) - start
      if snp == '-':
        seq[pos_in_seq] = ''
        delta_pos -= 1
      elif snp.startswith('I'):
        seq[pos_in_seq] = snp[1:]
        delta_pos += len(snp[2:])
      else:
        seq[pos_in_seq] = snp
      if not variant_subset or (contig, int(pos)) in variant_subset:
        if is_indel or snp != fasta.fetch(contig,int(pos),int(pos)+1):
          tr_snps[pos_in_seq+delta] = snp
      # if is_indel:
      delta += delta_pos
    out_seq = ''.join(seq)
    out_pos = ','.join(['%d:%s' % kv for kv in sorted(tr_snps.iteritems())])
    new_name = '%s-%d-%d-%s-%d' % (contig, start0, end0, fields[3], n)
    out.write('>%s\n%s\n' % (new_name, out_seq))
    out_variants.write('%s\t%s\n' % (new_name, out_pos))
    m += len(tr_snps)
    n += 1
