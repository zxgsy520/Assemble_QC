#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

#import taichi as ti
#ti.init()

LOG = logging.getLogger(__name__)

__version__ = "2.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []

#@ti.func
def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq = "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split("\n")
            yield seq[0], seq[1]
            seq = "%s\n" % line
        else:
            seq += line

    seq = seq.split("\n")
    if len(seq)==2:
        yield seq[0], seq[1]
    fp.close()


#@ti.func
def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip("@").split()[0]
            continue
        if line.startswith("@"):
            line = line.strip("@").split()[0]
            seq = seq.split("\n")

            yield seq[0], seq[1]
            seq = "%s\n" % line
        else:
            seq += "%s\n" % line

    if len(seq.split("\n"))==5:
        seq = seq.split("\n")
        yield seq[0], seq[1]

    fp.close()


#@ti.func
def split_scaffold(seq):    #分割序列

    tigseqs = []
    tigseq = ""
    gaplens = []
    n = 1
    tempn = 0

    for i in seq:
        n += 1
        if i == "N":
            if tigseq:
                end = n-2
                start = end - len(tigseq) + 1
                tigseqs.append(["%s-%s" % (start, end), tigseq])
                tigseq = ""
            tempn += 1
            continue
        if tempn:
            gaplens.append(tempn)
        start = n
        tempn = 0
        tigseq += i
    if tempn:
        gaplens.append(tempn)

    if tigseq:
        end = n - 1
        start = end - len(tigseq) + 1
        tigseqs.append(["%s-%s" % (start, end), tigseq])

    return gaplens, tigseqs


#@ti.func
def read_stat_seq(file, split=True):

    data = {}

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    data_scaf = {}
    tiglens = []
    gap = False
    scaflens = []

    for seq_id, seq in fh:
        seq = seq.upper()
        seqlen = len(seq)
        scaflens.append(seqlen)

        if "N" not in seq:
            data_scaf[seq_id] = [seqlen, []]
            tiglens.append(seqlen)
            if split:
                print(">%s\n%s" % (seq_id, seq))
        else:
            gap = True
            gaplens, tigseqs = split_scaffold(seq)
            data_scaf[seq_id] = [seqlen, gaplens]
            for site, seq in tigseqs:
                tiglens.append(len(seq))
                if split:
                    print(">%s:%s\n%s" % (seq_id, site, seq))

    return data_scaf, tiglens, scaflens, gap


def get_nx(lengths, x=0.5):

    addlen = 0
    total = sum(lengths)
    xl = 0

    for i in sorted(lengths, reverse=True):
        addlen += i
        if addlen >= total*x:
            xl = i
            break

    return xl
           

#@ti.func
def stat_contig(lengths):

    r = {"N50": [0, 0],
        "N60": [0, 0],
        "N70": [0, 0],
        "N80": [0, 0],
        "N90": [0, 0],
        "Longest": [0, 0],
        "Total": [0, 0],
        "Length>=1kb": [0, 0],
        "Length>=2kb": [0, 0],
        "Length>=5kb": [0, 0]
    }
    total = sum(lengths)
    maxlen = max(lengths)
    r["N50"][0] = get_nx(lengths, x=0.5)
    r["N60"][0] = get_nx(lengths, x=0.6)
    r["N70"][0] = get_nx(lengths, x=0.7)
    r["N80"][0] = get_nx(lengths, x=0.8)
    r["N90"][0] = get_nx(lengths, x=0.9)

    for i in lengths:
        if i >= 1000:
            r["Length>=1kb"][0] += i
            r["Length>=1kb"][1] += 1
        if i >= 2000:
            r["Length>=2kb"][0] += i
            r["Length>=2kb"][1] += 1
        if i >= 5000:
            r["Length>=5kb"][0] += i
            r["Length>=5kb"][1] += 1
        if i >= maxlen:
            r["Longest"][0] = i
            r["Longest"][1] += 1
        if  i >= r["N50"][0] :
            r["N50"][1] += 1
        if i >= r["N60"][0]:
            r["N60"][1] += 1
        if i >= r["N70"][0]:
            r["N70"][1] += 1
        if i >= r["N80"][0]:
            r["N80"][1] += 1
        if i >= r["N90"][0]:
            r["N90"][1] += 1

    r["Total"][0] = total
    r["Total"][1] = len(lengths)

    return r


#@ti.func
def stat_scaffold(data_scaf, scaflens):

    r = {"N50": [0, 0, 0, 0],
        "N60": [0, 0, 0, 0],
        "N70": [0, 0, 0, 0],
        "N80": [0, 0, 0, 0],
        "N90": [0, 0, 0, 0],
        "Longest": [0, 0, 0, 0],
        "Total": [0, 0, 0, 0],
        "Length>=1kb": [0, 0, 0, 0],
        "Length>=2kb": [0, 0, 0, 0],
        "Length>=5kb": [0, 0, 0, 0]
    }
    total = sum(scaflens)
    maxlen = max(scaflens)
    r["N50"][0] = get_nx(scaflens, x=0.5)
    r["N60"][0] = get_nx(scaflens, x=0.6)
    r["N70"][0] = get_nx(scaflens, x=0.7)
    r["N80"][0] = get_nx(scaflens, x=0.8)
    r["N90"][0] = get_nx(scaflens, x=0.9)
    r["Total"][0] = total
    r["Total"][1] = len(scaflens)

    for line in sorted(data_scaf.items(), key=lambda x:x[1], reverse=True):
        scaflen, gap = line[1]
        gapnum = len(gap)
        gaplen = sum(gap)

        r["Total"][3] += gapnum
        r["Total"][2] += gaplen

        if scaflen >= 1000:
            r["Length>=1kb"][0] += scaflen
            r["Length>=1kb"][1] += 1
            r["Length>=1kb"][3] += gapnum
            r["Length>=1kb"][2] += gaplen
        if scaflen >= 2000:
            r["Length>=2kb"][0] += scaflen
            r["Length>=2kb"][1] += 1
            r["Length>=2kb"][3] += gapnum
            r["Length>=2kb"][2] += gaplen
        if scaflen >= 5000:
            r["Length>=5kb"][0] += scaflen
            r["Length>=5kb"][1] += 1
            r["Length>=5kb"][3] += gapnum
            r["Length>=5kb"][2] += gaplen
        if scaflen >= maxlen:
            r["Longest"][0] = scaflen
            r["Longest"][1] += 1
            r["Longest"][2] = gaplen
        else:
            r["Longest"][3] += gapnum

        if scaflen >= r["N50"][0]:
            r["N50"][1] += 1
            r["N50"][3] += gapnum
            r["N50"][2] = gaplen
        if scaflen >= r["N60"][0]:
            r["N60"][1] += 1
            r["N60"][3] += gapnum
            r["N60"][2] = gaplen
        if scaflen >= r["N70"][0]:
            r["N70"][1] += 1
            r["N70"][3] += gapnum
            r["N70"][2] = gaplen
        if scaflen >= r["N80"][0]:
            r["N80"][1] += 1
            r["N80"][3] += gapnum
            r["N80"][2] = gaplen
        if scaflen >= r["N90"][0]:
            r["N90"][1] += 1
            r["N90"][3] += gapnum
            r["N90"][2] = gaplen

    return r


#@ti.func
def stat_genome(genome, outfile, split=True):

    data_scaf, tiglens, scaflens, gap = read_stat_seq(genome, split)

    r_tig = stat_contig(tiglens)
    fo = open(outfile, "w")
    if gap:
        fo.write("#StatType\tContigLength\tContigNumber\tScaffoldLength\tScaffoldNumber\tGapLength\tGapNumber\n")
        r_scaf = stat_scaffold(data_scaf, scaflens)
        for i in r_scaf:
            r_tig[i] += r_scaf[i]
    else:
        fo.write("#StatType\tContigLength\tContigNumber\n")

    for i in ["N50", "N60", "N70", "N80", "N90", "Longest", "Total", "Length>=1kb", "Length>=2kb", "Length>=5kb"]:
        temp = []
        for j in r_tig[i]:
            temp.append("{0:,}".format(j))

        fo.write("%s\t%s\n" % (i, "\t".join(temp)))

    fo.close()
    return 0


#@ti.func
def add_hlep_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help="Input genome sequence(fasta, fa, fast.gz)")
    parser.add_argument("-o", "--outfile", metavar="FILE", type=str, default="stat_genome.tsv",
        help="Input statistics output result file, default=stat_genome.tsv")
    parser.add_argument("-s", "--split", action="store_true", default=False,
        help=", default=False")

    return parser

#@ti.kernel
def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    stat_genome.py: statistical genome.
attention:
    stat_genome.py genome.fasta -o stat_genome.tsv
    stat_genome.py genome.fasta -o stat_genome.tsv --split >contig.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_genome(args.genome, args.outfile, args.split)


if __name__ == "__main__":

    main()
