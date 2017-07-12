#!/usr/bin/python

from thomtools import revComp
import subprocess
import argparse
import pysam
import sys
import re
import os

# parse command line
# parser = argparse.ArgumentParser(description = 'Format (and split) Stacks output fasta')
# parser.add_argument('-p','--stackpath', required = False,  default = './', help='path to Stacks output')
# parser.add_argument('-s','--samples',   required = True, help='TSV file of key/value sample pairs (sampleKey\tsampleVal\tpopulationVal)')
# parser.add_argument('-P','--sPlitBy', required = False, default = 'population', help='How should output fasta files be split?', choices = ['sample','population','none'])
# parser.add_argument('-shorten','--shorten_headers', required = False, default = 'no', help='Shorten fasta headers to only population_individual info?')
# parser.add_argument('-r','--prefix' , required = False, default = '' , help='prefix to add to output files')
# parser.add_argument('-o','--outPath', required = True,  default = './', help='Directory to output files')

# args=parser.parse_args()

# stacks  = args.stackpath
# prefix  = args.prefix
# sampsin = args.samples
# splitBy = args.sPlitBy
# shorten = args.shorten_headers
# outPath = args.outPath

# check if stacks fasta output exists
# fastain = ['batch_1.fa','batch_1.strict.fa']
# files  = os.listdir(stacks)
# if fastain[1] in files:
#     fastain = stacks+"/"+fastain[1]
#     print "using "+fastain+" as input."
# elif fastain[0] in files:
#     fastain = stacks+"/"+fastain[0]
#     print "using "+fastain+" as input."
# else:
#     sys.stderr.write("Required input file 'batch_1.fa' or 'batch_1.strict.fa' not present!")
 

#####################################
### DEFINE FUNCTIONS
#####################################

# PARSE FASTA FILE
def parseFA(faPATH):
    entries = {}
    faPATH = open(faPATH, 'r')
    while True:
        header = faPATH.readline()
        if not header: break
        header = header.strip('\n').lstrip(">")
        clocus = header
        seq    = faPATH.readline().strip('\n')
        entries[clocus] = seq
    faPATH.close()
    return entries

# MUTATE SNPS USING MD TAG IN SAM FILE
def mutateMD(seqArray, MDlist):
    nuc = 0
    mutated = seqArray
    numbers = re.compile("[0-9]")
    for mut in MDlist:
        if numbers.match(mut):
            # print mut+" exact matches"
            nuc += int(mut)
        elif "^" in mut: # SEEMS TO WORK FOR DELETIONS REL TO REF, SO FAR!!
            inslen = len(mut) - 1
        else:# NEED TO BE ABLE TO MUTATE ONE OR MORE NUCLEOTIDES
            snplist = list(mut)
            i = 0
            for s in snplist:
                # print "mutating "+seqArray[nuc+i]+" back to "+snplist[i]
                seqArray[nuc+i] = snplist[i]
                i += 1
            nuc += i
    return mutated

# ADD INSERTIONS RELATIVE TO REFERENCE USING CIGAR STRING
def mutateCigar(seqArray, ciglist):
    mutated = seqArray
    nuc = 0
    for cig in ciglist:
        nucs = int(cig[0:len(cig)-1])
        if "I" in cig: # I THINK THIS WILL DO IT
            for i in range(nucs):
                mutated[nuc+i] = "-"
            nuc += nucs
        elif 'D' not in cig:
            nuc += nucs
    return mutated

# PARSE BAM FILE OF ALIGNMENTS TO NINESPINE GENOME

# import sam file
samIN = "/share1/research/thom/seq/coal/chapter2/full/bbmap/stacks/rxstacks/phased_alaskanThreespines/aligned2nine/phasedtags.fa.mapPacBio.sam"
samIN = open(samIN, 'r')
sam   = {}

# DEFINE REGULAR EXPRESSIONS TO PARSE MD AND CIGAR TAGS
MDboundaries = "(?P<thing1>[A-Z]*)(?P<numbers>[0-9]+)(?P<thing2>[A-Z^]*)"
MDaddats     = "\g<thing1>@\g<numbers>@\g<thing2>"
Cboundaries  = "(?P<id>[0-9]+[=A-Z])"
Caddats      = "\g<id>@"
headermatch  = re.compile("^@[A-Z]+")
ntags = 0
print ""
for entry in samIN:
    if headermatch.match(entry):
        continue
    entry   = entry.strip('\n').split('\t')
    clocus  = entry[0]
    if entry[2] == "*":
        continue
    ntags+=1
    # if ntags < 20:
    #     print len(entry)
    if ntags % 1000 == 0:
        sys.stderr.write("Reading tag %s...\r"%(ntags))
    flag    = entry[1]
    subj    = entry[2]
    bpstart = entry[3]
    cigar   = entry[5]
    ciglist = re.sub(Cboundaries, Caddats, cigar)
    ciglist = ciglist[0:len(ciglist)-1].split("@")
    seq     = entry[9]
    MD      = entry[13]
    MDlist  = MD+":"
    MDlist  = re.sub(MDboundaries , MDaddats, MDlist)
    MDlist = MDlist.split("@")
    MDlen = len(MDlist)
    del MDlist[MDlen-1]
    del MDlist[0]
    mutated = list(seq)
    mutated = mutateMD(mutated, MDlist)
    mutated = mutateCigar(mutated, ciglist)
    if flag == '16':
        mutated = revComp("".join(mutated))
    else:
        mutated = "".join(mutated)
    # STORE AS DICTIONARY OF DICTIONARIES BECAUSE ASSOCIATIVE IS GOOOOOOOOOD!
    sam[clocus] = {'flag':flag, 'subj':subj, 'bp':bpstart,'cigar':cigar,'ciglist':ciglist,'seq':seq,'MD':MD, 'MDlist':MDlist, 'mutated':mutated}
samIN.close()

sys.stderr.write("  %s SAM alignments read.\n"%(ntags))

# FIND FASTA FILES

fastapath = "/share1/research/thom/seq/coal/chapter2/full/bbmap/stacks/rxstacks/phased_alaskanThreespines/fasta696/"
fastas    = os.listdir(fastapath)
fasOUT    = "/share1/research/thom/seq/coal/chapter2/full/bbmap/stacks/rxstacks/phased_alaskanThreespines/aligned2nine/fasta/"
sys.stderr.write("Finding fasta files and adding mutated ninespine consensus...\n")
nmatches = 0
for tag in sam:
    if tag+".fasta" in fastas:
        nmatches += 1
        alleles = parseFA(fastapath+tag+".fasta")
        faOUT = open(fasOUT+tag+".fasta", 'w')
        faOUT.write(">suqi\n")
        faOUT.write("%s\n"%(sam[tag]['mutated']))
        for allele in alleles:
            faOUT.write(">%s\n"%(allele))
            faOUT.write("%s\n"%(alleles[allele]))
        sys.stderr.write("  %s fasta files written...\r"%(nmatches))
sys.stderr.write("\n\n")







# for entry in sam:
#     print ""
#     print entry
#     print sam[entry]['cigar']
#     print sam[entry]['ciglist']
#     MD = sam[entry]['MD']
#     MDlist = sam[entry]['MDlist']
#     print MD+" ------------> "+" ".join(MDlist)
# print ""

