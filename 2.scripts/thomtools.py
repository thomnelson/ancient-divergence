#!/usr/bin/python
#
# thomtools.py by Thom Nelson
#
# Useful functions for performing coalescent simulations, generating sequence data,
#  and performing high-throughput analyses of sequence data. 
#

import random
import sys
import os
import re

###
### DEFINE FUNCTIONS TO RUN SIMULATIONS
###
def run_ms(nchroms, segsites, mu):
    # THIS FUNCTION IS SIMPLY A WRAPPER FOR RUNNING MS IN A PYTHON SCRIPT.
    # MS MUST BE IN THE PATH FOR THIS ONE TO WORK. 
    ms   = "ms {nsams} 1 -s {segs} -T -p 3 > ms.tmp"
    cmnd = ms.format(nsams = nchroms, segs = segsites)
    os.system(cmnd)
    ms_in = open("ms.tmp", "r")
    cmdline = ms_in.readline().strip('\n')
    rseed   = ms_in.readline().strip('\n')
    ms_in.readline().strip('\n'); ms_in.readline().strip('\n')
    newick  = ms_in.readline().strip('\n')
    nsegs   = ms_in.readline().strip('\n')
    positions = ms_in.readline().strip('\n').split(' '); positions = positions[1:-1]
    pos = []
    for position in positions:
        if "." in position:
            pos.append(position)
    haps = []
    for hap in ms_in:
        hap = hap.strip('\n')
        hap = list(hap)
        haps.append(hap)
    # Remove cutsite mutations if called for
    segsites = int(segsites)
    if mu == False:
        newpos  = []
        for i in range(len(pos)):
            test = int(float(pos[i]) * 1000)
            if test < 496 or test > 504:
                newpos.append(pos[i])
            else:
                segsites -= 1
                for j in range(len(haps)):
                    del haps[j][i]
        pos = newpos
        nsegs = "segsites: %s"%(segsites)
    for hap in range(len(haps)):
        haps[hap] = ''.join(haps[hap])
    ms = {'command':cmdline , 'random seed':rseed , 'tree':newick , \
          'segregating sites':str(nsegs) , 'positions':pos, 'haplotypes':haps}
    os.system('rm ms.tmp')
    return ms

def ms2halflotype(ms):
    # THIS FUNCTION TAKES HAPLOTYPES GENERATED WITH MS AND BREAKS THEM INTO TWO SETS OF ALLELES.
    #   THIS IS TO SIMULATED HAPLOTYPE DATA GENERATED WITH RAD, WHERE TWO ADJACENT HAPLOTYPES
    #   ARE SEQUENCED AND MUST BE PHASED.
    # THE OUTPUT IS A 'HALFOTYPES' CLASS, WHICH IS A DUMB NAME BUT CONTAINS ALL THE INFO NEEDED
    #   TO MATCH ALLELES TO INDIVIDUALS AND GENERATE INPUT FOR PHASE.
    unique  = []; halflos1 = []; halflos2 = []; indivs_bi = []; indivs_m = []; nalleles = 0
    hapLen  = len(ms['haplotypes'][1])
    sitesBefore = 0
    for pos in ms['positions']:
        pos = float(pos)
        if pos < 0.5:
            sitesBefore += 1
    
    cutsite = sitesBefore
    multi1  = {}
    translate1 = {}
    multi1n = 0
    multi2  = {}
    translate2 = {}
    multi2n = 0
    # populate lists of alleles present and dictionaries to relate biallelic halflotypes
    #   to their corresponding value in the 'multiallelic' genotype
    for allele in ms['haplotypes']:
        halflos = [allele[:cutsite],allele[cutsite:]]
        halflos1.append(halflos[0]); halflos2.append(halflos[1])
        if halflos[0] not in multi1:
            multi1n += 1
            multi1[halflos[0]] = multi1n
            translate1[multi1n] = halflos[0]
        if halflos[1] not in multi2:
            multi2n += 1
            multi2[halflos[1]] = multi2n
            translate2[multi2n] = halflos[1]
        if allele not in unique:
            unique.append(allele)
        nalleles += 1
    # populate indivs_* lists, which will be used to print out genotype inputs to PHASE
    i = 0
    while i < nalleles:
        geno1a_bi = halflos1[i] ; geno1b_bi = halflos1[i+1]
        geno2a_bi = halflos2[i] ; geno2b_bi = halflos2[i+1]
        geno1a_m  = multi1[geno1a_bi] ; geno1b_m = multi1[geno1b_bi]
        geno2a_m  = multi2[geno2a_bi] ; geno2b_m = multi2[geno2b_bi]

        indivs_bi.append([[geno1a_bi,geno2a_bi],[geno1b_bi, geno2b_bi]])
        indivs_m.append([[geno1a_m,geno2a_m],   [geno1b_m, geno2b_m]])
        i += 2

    # rescale positions so that RAD locus is 1 kb
    positions = ms['positions']
    for position in positions:
        position = float(position) * 1000
    # define 'unique' with true haplotypes and unique 'multiallelic' translations
    unique = {'true':unique, 'translate1':translate1, 'translate2':translate2}
    biallelic    = {'nindiv':nsam, 'nloci':str(hapLen),'positions':positions,'genotypes':indivs_bi}
    multiallelic = {'nindiv':nsam,'nalleles1':multi1n, 'nalleles2':multi2n , 'nloci':'2', 'positions':['0','1'],'genotypes':indivs_m}
    halflotypes  = {'unique':unique,'biallelic':biallelic,'multiallelic':multiallelic}
    return halflotypes


def generatePHASEinp(halflotypes, outfile_biallelic, outfile_multiallelic):
    # TAKE A 'HALFLOTYPES' OBJECT AND GENERATE A PHASE INPUT FILE
    bi    = halflotypes['biallelic']
    multi = halflotypes['multiallelic']
    
    # generate 'biallelic' output file
    bi_out = open(outfile_biallelic, 'w')
    indivs = bi['nindiv']+'\n' ; bi_out.write(indivs)
    nloci  = bi['nloci']+'\n'  ; bi_out.write(nloci)
    positions = bi['positions']
    for i in range(len(positions)):
        if float(positions[i]) > 0 and float(positions[i]) < 1:
            positions[i] = str(float(positions[i]) * 1000)
        else:
            positions[i] = str(positions[i])
    positions = 'P '+' '.join(positions)+'\n' ; bi_out.write(positions)
    segline = 'S'*int(bi['nloci'])+'\n' ; bi_out.write(segline)
    nindivs = 0
    for indiv in bi['genotypes']:
        nindivs += 1
        indivlabel = 'Sample_'+str(nindivs)+'\n'
        line1 = indiv[0][0]+indiv[0][1]+'\n'
        line2 = indiv[1][0]+indiv[1][1]+'\n'
        bi_out.write(indivlabel)
        bi_out.write(line1)
        bi_out.write(line2)
    bi_out.close()

    # generate 'multiallelic' output file
    m_out = open(outfile_multiallelic, 'w')
    indivs = multi['nindiv']+'\n' ; m_out.write(indivs)
    nloci  = multi['nloci']+'\n'  ; m_out.write(nloci)
    positions = 'P '+' '.join(multi['positions'])+'\n' ; m_out.write(positions)
    segline = 'M'*2+'\n' ; m_out.write(segline)
    nindivs = 0
    for indiv in multi['genotypes']:
        nindivs += 1
        indivlabel = 'Sample_'+str(nindivs)+'\n'
        line1 = str(indiv[0][0])+' '+str(indiv[0][1])+'\n'
        line2 = str(indiv[1][0])+' '+str(indiv[1][1])+'\n'
        m_out.write(indivlabel)
        m_out.write(line1)
        m_out.write(line2)
    m_out.close()

def runPHASE(inp,out,d, l, executable='phase'):
    # RUN PHASE, DUMP INTO A TEMPORARY FILE, AND COLLECT OUTPUT INTO A VARIABLE
    cmnd  = "{executable} -d{d} -l{l} -T {inp} {out} 1> PHASE.stdout 2> PHASE.stderr"
    PHASE = cmnd.format(executable = executable, d = d, l = l, inp = inp, out = out)
    os.system(PHASE)
    result = open(out, 'r')
    haplos = []
    for hap in result:
        hap = hap.strip('\n')
        haplos.append(hap)
    return haplos

###
### analyze PHASE output:
###  1) % correct best-guess haplos: compare to unique, ignore () or []
###  2) mean uncertainty: mean(len(haplo) - (segsites + (segsites - 1) / 2))
###  3) % correct @ indiv level: How many indivs have correctly called haplos?
###  4)
###

def sumPHASE(bi_out, m_out, halflotypes, cutsite):
    # CREATE A SUMMARY OF PHASE RESULTS FROM SIMULATED DATA
    resultS = [] ; resultM = []
    phaseS = open(bi_out,'r') ; phaseM = open(m_out,'r')
    # collect phaseS output
    alleles = 0 ; indivs  = 0
    guessesS = [] # collect n. phasing guesses per individual
    indiv   = [] # tmp holder for incoming alleles before addition of tuple to result
    for allele in phaseS:
        allele = allele.strip(' \n').split(' ')
        alleles += 1
        indiv.append(allele)
        if alleles % 2 == 0:
            resultS.append(indiv)
            guess = 0
            for site in allele:
                if '(' in site:
                    guess += 1
            guessesS.append(guess)
            indiv = []
            indivs += 1
    meanGuessesS = float(sum(guessesS)) / float(len(guessesS))
    # collect phaseM output
    alleles = 0 ; indivs  = 0; guessesM = []
    indiv   = [] # tmp holder for incoming alleles before addition of tuple to result
    for allele in phaseM:
        allele = allele.strip(' \n').split(' ')
        alleles += 1
        indiv.append(allele)
        if alleles % 2 == 0:
            resultM.append(indiv)
            guess = 0
            for site in allele:
                if '(' in site:
                    guess += 1
            guessesM.append(guess)
            indiv = []
            indivs += 1
    meanGuessesM = float(sum(guessesM)) / float(len(guessesM))
    # Compare phase output to original haplotypes
    trueS = halflotypes['biallelic'] ; trueM = halflotypes['multiallelic']
    
    # Because there are only two phasing 'options' per individual,
    #  confirm that first estimated haplotype matches either true haplotype
    #  if not, phasing is incorrect.
    
    # First, remove ambiguities '()' and '[]' from loci
    trueHaplosS = trueS['genotypes']
    for indiv in range(len(trueHaplosS)):
        for halflo in range(len(trueHaplosS[indiv])):
            full = ''.join(trueHaplosS[indiv][halflo])
            trueHaplosS[indiv][halflo] = full
    correctS = 0
    uniqueS = []
    for indiv in range(len(resultS)):
        for haplo in range(len(resultS[indiv])):
            for locus in range(len(resultS[indiv][haplo])):
                rmAmbig = re.sub('[\(\)\[\]]','',resultS[indiv][haplo][locus])
                resultS[indiv][haplo][locus] = rmAmbig
            haplotype = ''.join(resultS[indiv][haplo])
            if haplotype not in uniqueS:
                uniqueS.append(haplotype)
            resultS[indiv][haplo] = haplotype
        estimated = resultS[indiv]
        if estimated[0] in trueHaplosS[indiv]:
            correctS += 1

    trueHaplosM = trueM['genotypes']
    for indiv in range(len(trueHaplosM)):
        for haplo in range(len(trueHaplosM[indiv])):
            full = trueHaplosM[indiv][haplo]
            for half in range(len(full)):
                full[half] = str(full[half])
            trueHaplosM[indiv][haplo] = ''.join(full)

    correctM = 0
    uniqueM  = []
    for indiv in range(len(resultM)):
        for haplo in range(len(resultM[indiv])):
            full = ''.join(resultM[indiv][haplo])
            rmAmbig = re.sub('[\(\)\[\]]','',full)
            resultM[indiv][haplo] = rmAmbig
            if rmAmbig not in uniqueM:
                uniqueM.append(rmAmbig)
        estimated = resultM[indiv]
        if estimated[0] in trueHaplosM[indiv]:
            correctM += 1
    # collect results:
    #   chrom, coord, n. indiv, segsites, halflotype lengths, n. unique haplotypes, n unique halflotypes 1 and 2, 
    #   n est. haplotypes S and M, mean guessed phasings S and M, n. correct indivs S and M
    chrom = cutsite[0]
    coord = int(cutsite[1]) + 2
    lenHalf1 = 0
    for i in halflotypes['biallelic']['positions']:
        if float(i) < 500:
            lenHalf1 += 1
    lenHalf2 = int(segsites) - lenHalf1
    n_unique = len(halflotypes['unique'])
    n_unique1 = halflotypes['multiallelic']['nalleles1']
    n_unique2 = halflotypes['multiallelic']['nalleles2']
    results = [chrom, coord, int(nsam), int(segsites), lenHalf1, lenHalf2, n_unique, n_unique1, n_unique2, len(uniqueS), len(uniqueM),meanGuessesS,meanGuessesM,correctS,correctM]
    return results

def SAMheaders(nsam, SAMheader, out):
    for i in range(int(nsam)):
        fname   = out+"/Sample_"+str(i+1)+".sam"
        cmnd    = "cp {source} {dest}"
        SAMcopy = cmnd.format(source = SAMheader, dest = fname)
        os.system(SAMcopy)

def ms2sam(ms, coverage, nsam, fa, cutsite, out):
    ms = ms
    intpos = []
    # convert position strings into integers between 0-1000 and check if any hit cutsite
    for position in ms['positions']:
            position = int(float(position) * 1000)
            if position == 1000:
                position = 999
            intpos.append(position)
    # get sequence from stickleback genome
    chrom = cutsite[0] ; cutstart = int(cutsite[1]) ; cutend = int(cutsite[2]) ; cutlen = int(cutsite[3])
    midcoord = cutstart + (cutlen / 2)
    locusStart = midcoord - 500 ; locusEnd = midcoord + 499
    # Run samtools faidx to grab reference sequence
    cmnd = "samtools faidx {fasta} {chrom}:{locusStart}-{locusEnd} > tmp.fa"
    getAncestral = cmnd.format(fasta = fa, chrom = chrom, locusStart = str(locusStart), locusEnd = str(locusEnd))
    os.system(getAncestral)

    # Mutate ancestral sequence
    ancIN = open("tmp.fa", "r")
    ancIN.readline()
    ancestral = []
    for line in ancIN:
        line = line.strip('\n')
        ancestral.append(line)
    ancIN.close()
    os.system("rm tmp.fa")
    ancestral = ''.join(ancestral)
    ancestral = list(ancestral)
    ancSites = []
    derSites = []
    for i in intpos:
        anc = ancestral[i]
        ancSites.append(anc)
        derived = re.sub(anc, "" ,"ATCG") ; derived = random.sample(list(derived),1)[0]
        derSites.append(derived)
    mutations = [ancSites, derSites]

    # Generate mutated haplotypes and dump into SAM files
    haploCount  = 0
    sampleCount = 0
    sample      = []
    Samples     = []
    nsites      = len(intpos)
    for haplotype in ms['haplotypes']:
        allele = ancestral
        AncDir = list(haplotype)
        for site in range(len(AncDir)):
            AncDir[site] = int(AncDir[site])
        for i in range(nsites):
            position = intpos[i]
            state    = AncDir[i]
            base     = mutations[state][i]
            allele[position] = base
        allele = ''.join(allele)
        haploCount += 1
        
        if haploCount % 2 == 1:
            sample.append(allele)
        else:
            sample.append(allele)
            Samples.append(sample)
            sample = []
            sampleCount += 1
        rad1 = allele[0:502]
        rad2 = allele[498:]
    # Fake SAM fields
    image = 1000 ; x = 1000 ; y = 1000
    flag  = ['16','0'] ; mapq  = '38' ; cigar = '502M' ; rnext = '*' ; pnext = '0' ; tlen = '0' ; qual = 'I'*502
    for i in range(sampleCount):
        SAM = open(out+"/Sample_"+str(i+1)+".sam", "a")
        haplo1 = Samples[i][0] ; haplo2 = Samples[i][1]
        radA1  = haplo1[0:502] ; radA2  = haplo2[0:502]
        radB1  = haplo1[498:]  ; radB2  = haplo2[498:]
        halfsite1 = [radA1,radA2] ; halfsite2 = [radB1,radB2]
        halfsites = [halfsite1,halfsite2]
        coordA = locusStart
        coordB = locusStart + 498
        coords = [coordA,coordB]
        # Generate fake SAM fields        
        for locus in range(2):
            coord = str(coords[locus])
            strand = flag[locus]
            for allele in range(2):
                sequence = halfsites[locus][allele]
                for j in range(coverage):
                    flowcell = "1_%s_%s_%s_1"%(image,x,y)
                    entry = flowcell+"\t"+strand+"\t"+chrom+"\t"+coord+"\t"+mapq+"\t"+cigar+"\t"+rnext+"\t"+pnext+"\t"+tlen+"\t"+sequence+"\t"+qual+"\n"
                    SAM.write(entry)
                    image += 1 ; x += 1 ; y += 1
        SAM.close()


def readPopMap(popmap):
    # READ A POPMAP USED IN A STACKS ANALYSIS
    map      = {}
    popOrder = []
    popmapIN = open(popmap, 'r')
    for line in popmapIN:
        line = line.strip('\n').split('\t')
        pop = line[1]
        if pop in map:
            map[pop].append(line[0])
        else:
            map[pop] = [line[0]]
            popOrder.append(pop)
    popmapIN.close()
    return [map,popOrder]

def writeSubseq(file, fformat, subseqLen, taxa, outfile, popmap):
    # TRIM AN ALIGNMENT FILE TO A SPECIFIED LENGTH, OUTPUTTING A PHYLIP FILE
    popOrder = popmap[1]
    popmap   = popmap[0]
    npops    = len(popOrder)
    sampleN = len(taxa)
    aln = open(file, 'r')
    entries = {}
    if fformat == "phylip":
        header = aln.readline() ; header = header.strip('\n').split('\t')
        seqlen = header[-1]
        n      = header[-2]
        for entry in aln:
            entry = entry.strip('\n').split(' ')
            name = entry[0]
            seq = entry[-1]
            subseq = seq[0:subseqLen]
            entries[name] = subseq
    elif fformat == "fasta":
        l = 1
        name = ""
        for line in aln:
            if l % 2 == 1:
                name = line.strip('\n') ; name = name[1:]
            else:
                seq = line.strip('\n')
                subseq = seq[0:subseqLen]
                entries[name] = subseq
    aln.close()

    out = open(outfile, 'w')
    out.write(" %s\t%s\n"%(sampleN, subseqLen))
    for sample in taxa:
        seq = entries[sample]
        namelen = len(sample)
        gaplen  = 10 - namelen ; gaplen = " "*gaplen
        out.write("%s%s%s\n"%(sample,gaplen,seq))
    out.close()

def revComp(seq):
    # REVERSE COMPLEMENT A DNA SEQUENCE.
    #   REPLACE YOUR 'T's WITH 'U's LATER, RNASEEKERS!
    comp   = {'A':'T','T':'A','C':'G','G':'C','-':'-','N':'N','/':'/'}
    seq    = list(seq)
    seqlen = len(seq)
    seqnew = []
    for i in range(seqlen):
        pos  = seqlen - (i+1)
        base = seq[pos]
        seqnew.append(comp[base])
    return "".join(seqnew)
    
