#!/usr/bin/python
#
# Get per-individual and per-locus coverage for loci
#  that passed filtering in 'populations'
#

import argparse
import math
from gzip import GzipFile
import sys
import os

# parse command line options
parser = argparse.ArgumentParser(description='Extract Tmrca estimates from BEAST log files')
parser.add_argument('-haps','--haplotypes', required=False, default="./batch_1.haplotypes.tsv", help='haplotypes file from populations run')
parser.add_argument('-p','--popmap', required=True, help='Stacks population map for this run')
parser.add_argument('-d','--stackDIR', required=True, help='path to Stacks output files')

args=parser.parse_args()

hapfile = args.haplotypes
popmap  = args.popmap
stacks  = args.stackDIR
logfile = open("getStacksCoverage.log",'w')
logfile.write("script called using:\n  getStacksCoverage.py --haplotypes %s --popmap %s --stackDIR %s\n"%(hapfile,popmap,stacks))
# READ IN POPMAP TO ID INDIVIDUALS

samples = []
popmap = open(popmap, 'r')
for line in popmap:
    samples.append(line.strip("\n").split("\t")[0])
nsamples = len(samples)
sys.stderr.write("%s samples found...\n"%(nsamples))
logfile.write("%s samples found:\n"%(nsamples))
for sample in samples:
    logfile.write("%s "%(sample))
logfile.write("\n")
popmap.close()
# READ IN HAPLOTYPES FILE
# CREATE DICTIONARY TO STORE LOCUS AND COVERAGE FOR EACH SAMPLE (SET TO ZERO)
# THE INDEX OF THE SAMPLE IN THE 'SAMPLES' ARRAY WILL BE SAME AS THE INDEX OF
# THE SAMPLE IN THE VALUE OF LOCUS IN THE DICTIONARY
hapfile = open(hapfile, 'r')
locusCoverage = {}
hapfile.readline()
i = 0
sys.stderr.write("Reading haplotypes file...\n")
for tag in hapfile:
    i += 1
    tag = tag.split("\t")[0]
    locusCoverage[tag] = [0]*nsamples
hapfile.close()
sys.stderr.write("%s tags found in file.\n"%(i))
logfile.write("%s tags found in file.\n"%(i))

# LIST OUT STACKS DIRECTORY ENTRIES TO FIND MATCHES FILES AND FIGURE OUT IF GZIPPED

stackdir = os.listdir(stacks)
sys.stderr.write("Stacks output located in %s\n"%(stacks))
gzipped = 0
if samples[0]+".matches.tsv.gz" in stackdir:
    sys.stderr.write("Stacks output files appear to be gzipped.\n")
    logfile.write("Stacks output files appear to be gzipped.\n")
    gzipped = 1

# LOOP THROUGH SAMPLES.MATCHES FILES
# FOR EACH SAMPLE, READ IN MATCHES. CATALOG LOCUS IS [2], SAMPLE LOCUS IS [4], COV IS [6]
# ADD COVERAGE TO CORRECT INDEX FOR THE KEY == CATALOG LOCUS

for i in range(nsamples):
    sample  = samples[i]
    if gzipped == 1:
        matches = stacks+"/"+sample+".matches.tsv.gz"
        sys.stderr.write("Parsing "+sample+".matches.tsv.gz...\n")
        logfile.write("Parsing "+sample+".matches.tsv.gz...\n")
        matches = GzipFile(matches, 'r')
    else:
        matches = stacks+"/"+sample+".matches.tsv"
        sys.stderr.write("Parsing "+sample+".matches.tsv...\n")
        logfile.write("Parsing "+sample+".matches.tsv...\n")
        matches = open(matches, 'r')
    matches.readline()
    nmatches = 0
    for match in matches:
        match = match.strip('\n').split('\t')
        clocus = match[2]
        slocus = match[4]
        coverage = int(match[6])
        if clocus in locusCoverage:
            nmatches += 1
            locusCoverage[clocus][i] += coverage
            if nmatches % 1000 == 0:
                sys.stderr.write("Matches found (incl. alt. alleles):\t%s\r"%(nmatches))
    matches.close()
    sys.stderr.write("Matches found (incl. alt. alleles):\t%s\n"%(nmatches))
    logfile.write("Matches found (incl. alt. alleles):\t%s\n"%(nmatches))

# PRINT OUTPUT

outfile = open("stackCoverage.tsv",'w')
sys.stderr.write("Printing output to ./stackCoverage.tsv\n")
outfile.write("# script called using:\n#getStacksCoverage.py --haplotypes %s --popmap %s --stackDIR %s\n"%(hapfile,popmap,stacks))

samps = "\t".join(samples)
header = "locus"+"\t"+samps+"\n"
outfile.write(header)
for locus in locusCoverage:
    covs = locusCoverage[locus]
    for i in range(len(covs)):
        covs[i] = str(covs[i])
    outfile.write(locus+"\t"+"\t".join(covs)+"\n")
outfile.close()
logfile.close()
