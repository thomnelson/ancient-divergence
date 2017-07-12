#!/usr/bin/python
#
# stackPHASEr.py by Thom Nelson
#
# This script will take Stacks output files (after running populations)
#   and use PHASE to attempt to phase adjacent RAD loci. Stacks must be run
#   with either pstacks or through ref_map.pl to identify adjacent tags. It is
#   definitely possible to phase tags based on LD calculations or that were
#   aligned post-hoc, but this script won't do that. I'd be more than happy to
#   write one for those cases later.
#

from thomtools import *
import argparse
import random
import sys
import os
import re

#############################
# parse command line options
#############################
parser = argparse.ArgumentParser(description='Use PHASE to phase adjacent RAD loci into full RAplotypes')
# required options
parser.add_argument('-p','--path', required=True, help='Path to Stacks output')
parser.add_argument('-b','--batch', required=False,default = '1', help='Batch ID (default = 1)')
parser.add_argument('-m','--popmap', required=False,default = '', help='Population map from Stacks. Used to shorten/specify sample names')
parser.add_argument('-o','--out', required=False,default = '||', help='Path to write output files (default = /stack/path/PHASEd')

args=parser.parse_args()
stackpath = args.path
batch     = args.batch
popmap    = args.popmap
outpath   = args.out
if outpath == '||':
    outpath = stackpath+'/PHASEd/'
else:
    outpath = outpath+"/"


#############################
# Check if files and paths exist. Make directories if need be
#############################

if os.path.exists(outpath) == False:
    os.system("mkdir "+outpath)

stacks_output = os.listdir(stackpath)

needed = ['batch_1.catalog.tags.tsv','batch_1.catalog.snps.tsv','batch_1.haplotypes.tsv']
notThere = []
crash = 0
for f in needed:
    if f not in stacks_output:
        crash = 1
        notThere.append(f)
if crash == 1:
    sys.stderr.write('ERROR: The following files are not present, or may be compressed:\n')
    for f in notThere:
        sys.stderr.write(f+'\n')
    sys.exit()

tags        = stackpath+'/'+needed[0]
snps        = stackpath+'/'+needed[1]
haplotypes  = stackpath+'/'+needed[2]
log         = open(outpath+"/stackPHASEr.log", 'w')
pairlog     = open(outpath+"/stackPHASEr.pairs", 'w')

log.write("call: stackPHASEr.py -p "+stackpath+" -b "+batch+" -m "+popmap+" -o "+outpath+"\n")
log.write("Input files:\n"+tags+"\n"+haplotypes+"\n")

#############################
# Create function to convert alleles (segsites only) to
#  full-length haplotypes for export
#############################

def alleles2fasta(consensus, columns, alleles, sampleNames, faOUT):
    consensus = list(consensus)
    faOUT = open(faOUT,'w')
    for s in range(len(alleles)):
        sampName = sampleNames[s]
        haps = alleles[s] ; haps = haps.split('/')
        if len(haps) == 1:
            haps = ''.join(haps)
            haps = [haps,haps]
        for i in range(len(haps)):
            if "-" not in haps[i]:
                header = ">"+sampName+"_"+str(i)+"\n"
                sites = haps[i]
                hap   = consensus
                for j in range(len(sites)):
                    col = columns[j]
                    base = sites[j]
                    hap[col] = base
                haps[i] = "".join(hap)
                seq = haps[i]
                faOUT.write(header+seq+"\n")

def phaseMulti2alleles(phased, halflotypes, locus = '0'):
    phased = phased
    trans1 = halflotypes['unique']['translate1']
    trans2 = halflotypes['unique']['translate2']
    consensus = halflotypes['consensus']
    positions = halflotypes['biallelic']['positions']
    for p in range(len(positions)):
        positions[p] = int(positions[p])
    nsamps = 0
    nhaps  = 0
    sample  = []
    samples = []
    #####################################
    # Read through phased haplotypes and
    # translate back to nucleotides
    #####################################
    for haplo in range(len(phased)):
        hap = re.sub("[\(\)\[\]]", "", phased[haplo])
        hap = hap.split(' ')
        if int(hap[0]) in trans1:
            hap[0] = trans1[int(hap[0])]
        else:
            hap[0] = '-'
        if int(hap[1]) in trans2:
            hap[1] = trans2[int(hap[1])]
        else:
            hap[1] = '-'
        hap = ''.join([hap[0],hap[1]])
        nhaps += 1
        if nhaps % 2 == 1:
            sample.append(hap)
        else:
            sample.append(hap)
            sample = "/".join(sample)
            samples.append(sample)
            sample = []
            nsamps += 1
    s = 1
    missingdata = []
    locusOut = [consensus, positions, samples]
    for sample in samples:
        haplo1 = sample[0]
        haplo2 = sample[1]
        if '-' in haplo1:
            missing = 'Locus %s, allele 1, for sample %s is unknown and will not be exported.\n'%(locus,s)
            missingdata.append(missing)
        if '-' in haplo2:
            missing = 'Locus %s, allele 2, for sample %s is unknown and will not be exported.\n'%(locus,s)
            missingdata.append(missing)
        s += 1
    locusOut = [locusOut, missingdata]
    return locusOut

######################            
# Read in catalog SNPs file to extract positions of segregating sites w/i haplotype
######################

sys.stderr.write("Reading SNP positions in from %s...\n"%(needed[1]))
columns = {}
snps = open(snps, 'r')
snps.readline()
nsnps = 0
for snp in snps:
    nsnps += 1
    snp    = snp.strip("\n").split("\t")
    locus  = snp[2] ; col = snp[3]
    if locus in columns:
        columns[locus].append(int(col))
    else:
        columns[locus] = [int(col)]
snps.close()
sys.stderr.write("SNPs found in the catalog:\t%s\n"%(nsnps))
log.write("SNPs found in the catalog:\t%s\n"%(nsnps))
######################            
# Read in catalog tags file to generate locusInfo
######################

sys.stderr.write('Reading %s...\n'%(needed[0]))
lociPlus  = {}
lociMinus = []
tagsIN = open(tags,'r')
tagsIN.readline()

print "reading tags"

for tag in tagsIN:
    tag = tag.strip('\n').split('\t')
    clocus = tag[2]; group = tag[3] ; bp = tag[4] ; strand = tag[5] ; consensus = tag[9]
    cols = []
    if clocus in columns:
        cols   = columns[clocus]
    locus  = [clocus, group, bp, strand, consensus, cols]
    if strand == "+":
        key = group+bp+strand
        val = {'locus':clocus,'consensus':consensus,'columns':cols}
        lociPlus[key] = val
    else:
        lociMinus.append(locus)
tagsIN.close()

#####################
# If locus is on minus strand, reverse complement seq, clip off leading 'TGCA',
#  change 'column' assignments, and check to see if its pair is present in lociPlus 
#####################

pairs = {}
found = 0
print "looking at loci on minus strand"
        
for i in range(len(lociMinus)):
    locus  = lociMinus[i]
    clocus = locus[0] ; group = locus[1] ; bp = locus[2] ; consensus = locus[4] ; cols = locus[5]
    overhang = 4
    consensus = consensus[overhang:]
    consensus = revComp(consensus)
    seqlen = len(consensus)
    for j in range(len(cols)):
        cols[j] = seqlen - cols[j]
    lociMinus[i][5] = cols
    # find expected offset of locus pair to search for matching bp position in pair
    offset = int(bp) - 3 ; offset = str(offset)
    keycheck = group+offset+'+'
    if keycheck in lociPlus:
        pair = lociPlus[keycheck]
        pairtext = clocus+"\t"+pair['locus']+"\n"
        pairlog.write(pairtext)
        pairs[clocus] = {'group':group, 'bp':bp, 'pair':pair ,'consensus':consensus, 'columns':cols}
        found += 1
        sys.stderr.write("Locus pairs found:\t%s\r"%(found))

log.write("Locus pairs found:\t%s\n"%(found))
print "Locus pairs found:\t%s"%(found)
print "Grabbing filtered genotypes from %s..."%(needed[2])

################################################
# Haplotype file will be more stringent than catalog files, since
# this is produced by populations. This, then, needs some checking/cross-referencing.
#
# Additionally, grab sample names from header and shorten if necessary/desired
################################################

hapsIN = open(haplotypes,'r')
hapHeader = hapsIN.readline()
hapHeader = hapHeader.strip('\n').split('\t') ; hapHeader = hapHeader[2:]

sampleNames = hapHeader

if len(popmap) > 0:
    popmap = open(popmap,'r')
    i = 0
    s = 0
    last = ""
    for line in popmap:
        line = line.strip('\n').split('\t')
        if i == 0:
            last = line[1]
        current = line[1]
        if current == last:
            s += 1
        else:
            s = 1
        n = "{:0>2d}".format(s)
        sampleNames[i] = line[1]+"_"+str(n)
        last = current
        i += 1

halfs = {}
i = 1
for locus in hapsIN:
    locus = locus.strip('\n').split('\t')
    clocus = locus[0] ; genotyped = locus[1] ; genotypes = locus[2:]
    halfs[clocus] = genotypes
    i += 1
hapsIN.close()
print "RAD loci present after population filtering:\t%s"%(i)
log.write("RAD loci present after population filtering:\t%s\n"%(i))


passedFilter = 0  # increment if both locus and pair showed up in 'haplotypes' file
nopass = []
noPair = 0
merge  = 0

################################################
# assign genotypes to existing locus pairs
# If one locus in pair doesn't exist or is 'consensus',
# print remaining locus directly to a fasta file.
################################################

for locus in pairs:
    noSNP      = [0,0] # add 1 each time 'consensus' shows up. write as fa if sum == 1
    nopair     = 0
    unphased   = {}
    Pass       = True
    pair       = pairs[locus]['pair']
    pairlocus  = pair['locus']
    genosMinus = []
    genosPlus  = []
    consensusMinus = pairs[locus]['consensus']
    consensusPlus  = ''
    columnsMinus   = pairs[locus]['columns']
    columnsPlus    = []
    # locus is (-)strand in pair
    if locus in halfs:
        # Locus has called genotypes (even if no SNPs)
        # Evaluate checks to decide how to proceed
        ###############################
        genosMinus = halfs[locus]
        # Check if 'minus' halfsite contains SNPs
        if 'consensus' in genosMinus:
            noSNP[0] = 1 ; Pass = False
        else:
            pairs[locus]['genos(-)'] = genosMinus
        # Check if paired locus has called genotypes and SNPs
        if pairlocus in halfs:
            genosPlus = halfs[pairlocus]
            consensusPlus = pair['consensus']
            columnsPlus   = pairs[locus]['pair']['columns']
            if 'consensus' in genosPlus:
                noSNP[1] = 1 ; merge += 1 ; Pass = False
            else:
                pairs[locus]['genos(+)'] = genosPlus
                passedFilter += 1 ; merge += 1
                # Loci that 'passedFilter' will be sent to PHASE
        else:
            nopair += 1 ; noPair += 1 ; Pass = False
    # If no genotypes called for minus-strand locus, try +-strand locus
    elif pairlocus in halfs:
        consensusPlus = pairs[locus]['pair']['consensus']
        columnsPlus   = pairs[locus]['pair']['columns']
        nopair += 1 ; noPair += 1 ; noSNP[0] = 1
        genosPlus = halfs[pairlocus]
        if 'consensus' in genosPlus:
            noSNP[1] = 1
        Pass = False
    else:
        Pass = False
    if Pass == False:
        nopass.append(locus)

    # Use noSNP values to sort output
    # noSNP  = [1,1]          : Both halfsites are monomorphic. Makes for easy exporting!
    # noSNP  = [0,1] or [1,0] : one halfsite in pair is missing or monomorphic
    # nopair > 0              : only one halfsite was present or passed Stacks filtering
    if sum(noSNP) == 1 and nopair == 0:
        stdout = "Loci {l1} and {l2} merged without phasing. Locus {l3} contained no SNPs\n"
        # One halfsite contains no SNPs. Mutate other halfsite and join to form full locus
        consensus = consensusMinus + consensusPlus
        if noSNP[0] == 1:
            shift = len(consensusMinus)
            for i in range(len(columnsPlus)):
                columnsPlus[i] = columnsPlus[i] + shift
            faOUT = outpath+pairlocus+".fasta"
            alleles2fasta(consensus = consensus, columns = columnsPlus, alleles = genosPlus, sampleNames = sampleNames, faOUT = faOUT)
            report = stdout.format(l1 = locus, l2 = pairlocus, l3 = locus)
            log.write(report)
        else:
            faOUT = outpath+locus+".fasta"
            alleles2fasta(consensus = consensus, columns = columnsMinus, alleles = genosMinus, sampleNames = sampleNames, faOUT = faOUT)
            report = stdout.format(l1 = locus, l2 = pairlocus, l3 = locus)
            log.write(report)
    if sum(noSNP) == 1 and nopair > 0:
        stdout = "No genotypes found for locus {l1}. Locus {l2} exported by its lonesome.\n"
        if noSNP[0] == 1:
            faOUT = outpath+pairlocus+".fasta"
            alleles2fasta(consensus = consensusPlus, columns = columnsPlus, alleles = genosPlus, sampleNames = sampleNames, faOUT = faOUT)
            report = stdout.format(l1 = locus, l2 = pairlocus)
            log.write(report)
        else:
            faOUT = outpath+locus+".fasta"
            alleles2fasta(consensus = consensusMinus, columns = columnsMinus, alleles = genosMinus, sampleNames = sampleNames, faOUT = faOUT)
            report = stdout.format(l1 = pairlocus, l2 = locus)
            log.write(report)
    #
    # Need an output when no polymorphisms detected in either/both half sites
    #
    if sum(noSNP) == 2:
        if nopair == 0:
            # simply merge consensus seqs
            consensus = consensusMinus+consensusPlus
            faOUT = outpath+locus+".fa"
            alleles2fasta(consensus=consensus,columns=[],alleles=[], sampleNames = sampleNames,faOUT=faOUT)
            stdout = "Both halfsites found for loci {l1} and {l2}. No polymorphisms detected in either.\n"
            report = stdout.format(l1 = locus, l2 = pairlocus)
            log.write(report)
        elif 'consensus' in genosPlus:
            consensus = consensusPlus
            faOUT = outpath+pairlocus+".fa"
            alleles2fasta(consensus = consensus,columns=[],alleles=[], sampleNames = sampleNames,faOUT=faOUT)
            stdout = "Locus {l} has no pair and no polymorphisms. How sad...\n"
            report = stdout.format(l = pairlocus)
            log.write(report)
        else:
            consensus = consensusMinus
            faOUT = outpath+locus+".fa"
            alleles2fasta(consensus = consensus,columns=[],alleles=[], sampleNames = sampleNames,faOUT=faOUT)
            stdout = "Locus {l} has no pair and no polymorphisms. How sad...\n"
            report = stdout.format(l = locus)
            log.write(report)
print "Clearing merged loci from memory..."
for locus in nopass:
    del pairs[locus]
    
log.write("Valid loci without pairs:\t%s\n"%(noPair))
log.write("Valid pairs merged without phasing:\t%s\n"%(merge))
log.write("Valid pairs to be phased:\t%s\n"%(passedFilter))

print "Valid loci without pairs:\t%s"%(noPair)
print "Valid pairs merged without phasing:\t%s"%(merge)
print "Valid pairs to be phased:\t%s"%(passedFilter)

# We now have a collection of all loci that passed filtering and should be phased
# pairs = {group, bp, pair, genos-, genos+}
# Generate halflotype object?
print "Generating halflotype objects..."
halflotypesAll = {}
for locus in pairs:
    genosMinus  = pairs[locus]['genos(-)']
    nindivMinus = len(genosMinus)
    genosPlus   = pairs[locus]['genos(+)']
    nindivPlus  = len(genosPlus)
    # Set up consensus sequence and columns (e.g. 'positions') for whole locus
    consensus    = pairs[locus]['consensus'] + pairs[locus]['pair']['consensus']
    shift        = len(pairs[locus]['consensus'])
    columnsMinus = pairs[locus]['columns']
    columnsPlus  = pairs[locus]['pair']['columns']
    positions = []
    for geno in genosMinus:
        if '-' not in geno:
            ncols = len(geno)
            if '/' in geno:
                gen = geno.split('/')
                ncols = len(gen[0])
            for i in range(ncols):
                positions.append(columnsMinus[i])
            break
    for geno in genosPlus:
        if '-' not in geno:
            ncols = len(geno)
            if '/' in geno:
                gen = geno.split('/')
                ncols = len(gen[0])
            for i in range(ncols):
                positions.append(str(columnsPlus[i] + shift))
            break
        
    missingMinus= [] ; missingPlus = []
    lenMinus    = 0  ; lenPlus     = 0
    # Split alleles and identify missing data
    for i in range(nindivMinus):
        genosMinus[i] = revComp(genosMinus[i])
        if '/' in genosMinus[i]:
            genosMinus[i] = genosMinus[i].split('/')
            lenMinus = len(genosMinus[i][0])
        elif genosMinus[i] == '-':
            missingMinus.append(i)
            genosMinus[i] = ['-','-']
        else:
            genosMinus[i] = [genosMinus[i],genosMinus[i]]
            lenMinus = len(genosMinus[i][0])

    for i in range(nindivPlus):
        if '/' in genosPlus[i]:
            genosPlus[i] = genosPlus[i].split('/')
            lenPlus = len(genosPlus[i][0])
        elif genosPlus[i] == '-':
            missingPlus.append(i)
            genosPlus[i] = ['-','-']
        else:
            genosPlus[i] = [genosPlus[i],genosPlus[i]]
            lenPlus = len(genosPlus[i][0])
    ##################################
    # Start creating halflotype object
    ##################################
    nindiv = nindivMinus
    nloci = lenMinus + lenPlus
    positionsB = positions
    positionsM = ['0','1']
    genotypesB = []
    genotypesM = []
    multiMinus = {}
    multiPlus  = {}
    nuniqMinus = 0
    nuniqPlus  = 0
    translateMinus = {} # translate multiallelic genotypes back to SNPs
    translatePlus  = {}
    for i in range(nindiv):
        minus = genosMinus[i]
        plus  = genosPlus[i]
        minusM= []          # hold multiallelic genotypes
        plusM = []
        for j in range(2):
            if minus[j] in multiMinus:
                minusM.append(multiMinus[minus[j]])
            elif minus[j] == '-':
                minus[j] = '?'*lenMinus
                minusM.append('-1')
            else:
                nuniqMinus += 1
                minusM.append(nuniqMinus)
                multiMinus[minus[j]] = nuniqMinus
                translateMinus[nuniqMinus] = minus[j]
            if plus[j] in multiPlus:
                plusM.append(multiPlus[plus[j]])
            elif plus[j] == '-':
                plus[j] = '?'*lenPlus
                plusM.append('-1')
            else:
                nuniqPlus += 1
                plusM.append(nuniqPlus)
                multiPlus[plus[j]] = nuniqPlus
                translatePlus[nuniqPlus] = plus[j]
        exportB = [[minus[0],plus[0]],[minus[1],plus[1]]]
        exportM = [[minusM[0],plusM[0]],[minusM[1],plusM[1]]]
        genotypesB.append(exportB)
        genotypesM.append(exportM)
    # Format data to create a 'halflotypes' object
    biallelic             = {'nindiv': str(nindiv),'nloci': str(nloci),'positions': positionsB,'genotypes': genotypesB}
    multiallelic          = {'nindiv': str(nindiv),'nloci': '2','positions': positionsM,'genotypes': genotypesM}
    unique                = {'translate1':translateMinus,'translate2':translatePlus}
    consensus             = consensus
    halflotypesAll[locus] = {'biallelic':biallelic,'multiallelic':multiallelic,'unique':unique, 'consensus':consensus}



loci = 0
for halflotype in halflotypesAll:
    halflotypes = halflotypesAll[halflotype]
    sys.stderr.write("                                                                              \r")
    sys.stderr.write("Curently phasing catalog locus %s. Total loci phased:\t%s\r"%(halflotype,loci))

    outB = outpath+halflotype+"b.inp"
    outM = outpath+halflotype+"m.inp"
    generatePHASEinp(halflotypes = halflotypes, outfile_biallelic = outB, outfile_multiallelic = outM)
    phased = runPHASE(inp = outM, out = 'tmp.out', d = 1, l = 2)
    alleles = phaseMulti2alleles(phased = phased, halflotypes = halflotypes, locus = halflotype)
    errors    = alleles[1]
    consensus = alleles[0][0]
    positions = alleles[0][1]
    alleles   = alleles[0][2]
    for missing in errors:
        log.write(missing)
    alleles2fasta(consensus = consensus, columns = positions, alleles = alleles, sampleNames = sampleNames, faOUT = outpath+halflotype+".fasta")
    loci += 1
    
os.system('rm tmp.out*')
sys.stderr.write('\n')
