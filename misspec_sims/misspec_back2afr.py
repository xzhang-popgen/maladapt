#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 17:09:00 2022

@author: xinjunzhang
"""



import msprime, pyslim, os, random, itertools,argparse
import numpy as np

parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-s', '--sim', action="store", dest="sim_id",
                        help="which simulation batch, default: 1",
                        default=1, type=int)
parser.add_argument('-g', '--gene', action="store", dest="gene_id",
                        help="which simulation batch, default: 1; range 1-26",
                        default=1, type=int)
parser.add_argument('-d', '--dominance', action="store", dest="dominance_id",
                        help="which simulation batch, default: 1; range 1-26",
                        default=1, type=int)
args = parser.parse_args()

batch = str(args.sim_id) #now sim id is going to refer to a 
whichgene = int(args.gene_id) 
whichh = int(args.dominance_id) 

######################################################################
def getSelectedAlleleSelcoeff(file_name): # returns true if m2 was fixed in p1 before admixture
    infile = open(file_name)
    m2_in_p4 = 0 
    for line in infile:
        if line[0:5]=='#OUT:':
            fields = line.split()
            selcoeff = fields[7]
            if (int(fields[1]) == 8900) & (fields[3] == "p4"):
            	m2_in_p4 = 1
            #break
    infile.close()
    return selcoeff,m2_in_p4

def simplify_tree(filename, sample_sizes, archaic_sample_time):
    '''
    this method subsamples a tree sequence to a given set of
    sample sizes. archaic sample is set to be 1 individual in SLiM
    '''
    #load tree sequence
    ts = pyslim.load(filename)
    indivs = [x for x in ts.individuals()]
    #get nodes (chromosomes) alive today
    nodes_today_p1 = [x.id for x in ts.nodes() if ((x.time == 0.0) and
                      (x.population == 1))]
    nodes_today_p4 = [x.id for x in ts.nodes() if ((x.time == 0.0) and
                      (x.population == 4))]
    #get a list of the individuals alive today
    indivs_today_p1 = [x.id for x in indivs if x.nodes[0] in
                       nodes_today_p1]
    indivs_today_p4 = [x.id for x in indivs if x.nodes[0] in
                       nodes_today_p4]
    #subsample individuals to sample sizes
    indivs_sample_p1 = random.sample(indivs_today_p1, sampsize[0])
    indivs_sample_p4 = random.sample(indivs_today_p4, sampsize[1])
    #get their nodes
    nodes_sample_p1 = [x.nodes for x in indivs if x.id in
                       indivs_sample_p1]
    nodes_sample_p4 = [x.nodes for x in indivs if x.id in
                       indivs_sample_p4]
    nodes_sample_p1 = list(itertools.chain.from_iterable(
                           nodes_sample_p1))
    nodes_sample_p4 = list(itertools.chain.from_iterable(
                           nodes_sample_p4))
    #get archaic samples
    archaic = [x.id for x in ts.nodes() if x.time ==
               archaic_sample_time]
    indivs_archaic = [x.id for x in indivs if x.nodes[0] in archaic]
    nodes_sample_archaic = [x.nodes for x in indivs if x.id in
                            indivs_archaic]
    nodes_sample_archaic = list(itertools.chain.from_iterable(
                                nodes_sample_archaic))
    #subsample while retaining admixture recorded nodes
    samp = nodes_sample_archaic + nodes_sample_p1 + nodes_sample_p4
    ts_sample = ts.simplify(samples=samp, filter_populations=False)
    return ts_sample

def ns_vcf(ts, vcfoutpath):
    with open(vcfoutpath, "w") as vcf_file:
        ts.write_vcf(vcf_file, ploidy=2)

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1
            mut = site.mutations[0]
            keep_mutation = True
            for i in range(len(start)):
                left = start[i]
                right = end[i]
                assert(left < right)
                if i > 0:
                    assert(end[i - 1] <= left)
                if left <= site.position < right:
                    keep_mutation = (random.uniform(0, 1) > proportion)
            if keep_mutation:
                site_id = tables.sites.add_row(
                    position=site.position,
                    ancestral_state=site.ancestral_state)
                tables.mutations.add_row(
                    site=site_id, node=mut.node, derived_state=mut.derived_state)
    return tables.tree_sequence()

def noncoding_vcf(ts, vcfoutpath, mu, start, end, proportion=1.):
    ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=float(1.8e-8),
         keep=False))
    proportion = 1.
    ts = remove_mutations(ts, start, end, proportion)
    with open(vcfoutpath, "w") as vcf_file:
        ts.write_vcf(vcf_file, ploidy=2)

def syn_vcf(ts, vcfoutpath, mu, start, end):
    mu = 1/3.31*mu
    ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=float(mu),
         keep=False))
    proportion = 1.0
    first = ts.first().interval[0]
    last = ts.last().interval[1]
    new_start = [first] + [x for x in end[0:(len(end))]]
    new_end = [x for x in start[0:(len(start))]] + [last]
    ts = remove_mutations(ts, new_start, new_end, proportion)
    with open(vcfoutpath, "w") as vcf_file:
        ts.write_vcf(vcf_file, ploidy=2)


##################################################

def calc_p1ancestry (treepath, admpop, popsize,t_sinceadm):
    #treepath = "test.trees"
    ts = pyslim.load(treepath)
    any_ancestry = ancestry_p_varies(ts,admpop,popsize,t_sinceadm)
    meanp1 = sum(any_ancestry)/len(any_ancestry)
    return meanp1

def ancestry_p_varies(ts,pop,nsize,duration): #pop=source pop
    #nsize = 4108
    n=nsize*2
    mixtime=duration+1
    #times = [x.time for x in ts.nodes() if (x.population == 2)]
    #ids = [x.id for x in ts.nodes() if (x.population == 2)]
    #times = list(set(times))
    #times.sort()
    p = [x.id for x in ts.nodes() if ((x.population == int(pop)) and (x.time == mixtime))] #source pop
    today = [x.id for x in ts.nodes() if ((x.population == 4) and (x.time == 0))] #assuming p3 is recipient
    tree_p = [sum([t.num_tracked_samples(u) for u in p])/n
               for t in ts.trees(tracked_samples=today, sample_counts=True)]      
    starts=[]
    ends=[]
    for x in ts.trees():
        starts.append(x.interval[0])
        ends.append(x.interval[1])
    #sum(tree_p)/len(tree_p)
    return tree_p


def vcf2genos (vcfpath,p1size,p2size,p3size):   #sizes are int
    #vcfpath = "test_combined.vcf"
    #p1size,p2size,p3size = 2,100,100
    vcf = open(vcfpath)
    start = 0
    startg = 0
    pos = []
    geno = []
    muttype = []
    for count,line in enumerate(vcf):
        if line[0:6]=="#CHROM":
            start +=1
            startg = count
        if (start == 1) & (startg!=0) & (count>startg):
            fields = line.split()
            if ("1|2" not in fields) and ("0|2" not in fields) and ("2|0" not in fields) and ("2|1" not in fields) and ("2|2" not in fields):
                pos.append(int(fields[1]))
                muttype.append(fields[7][3:])
                geno.append(fields[9:])
    vcf.close()
    len_genome = int(pos[len(pos)-1])-int(pos[0])    
    geno1 = []
    geno2 = []
    geno3 = []
    for n in range(len(geno)):
        geno1.append([geno[n][x] for x in range(p1size)])
        geno2.append([geno[n][x] for x in range(p1size,p1size+p2size)])
        geno3.append([geno[n][x] for x in range(p1size+p2size,p1size+p2size+p3size)])
    return geno1, geno2, geno3, len_genome, muttype, pos

def geno2hap (geno):
    #geno=geno1    
    num_pos = len(geno)
    num_ind = len(geno[1])    
    geno_ref = np.empty([num_pos,num_ind])
    geno_alt = np.empty([num_pos,num_ind])
    for n in range(num_pos):
        geno_ref[n] = [int(i.split('|', 1)[0]) for i in geno[n]]
        geno_alt[n] = [int(i.split('|', 1)[1]) for i in geno[n]]
    geno_ref = np.transpose(geno_ref)
    geno_alt = np.transpose(geno_alt)
    hapmat = np.append(geno_ref,geno_alt,axis=0)
    return hapmat

def calc_freq (pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq/ float(pop_hap.shape[0])
    return popfreq


def vSumFunc(other_hap, currentArchi):
    current_hap = np.array([p1_hap[currentArchi,]])
    div = np.zeros(other_hap.shape)
    ones = np.ones((other_hap.shape[0],1))
    current_hap = current_hap
    current_hap_extended = np.dot(ones, current_hap)
        #print(current_hap_extended)
        #computes vectorized logical xor between current_hap and the hap being analyzed
    div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
        #reduces the div on each row with sum to count the number of different alleles
    return np.add.reduce(div, 1)   


def calc_stats (p1_hap,p2_hap,p3_hap,len_genome):
    #D
    p1_freq = calc_freq (p1_hap)
    p2_freq = calc_freq (p2_hap)
    p3_freq = calc_freq (p3_hap)
    abbavec = (1.0 - p2_freq)*p3_freq*p1_freq
    babavec = p2_freq*(1.0 - p3_freq)*p1_freq
    abba = np.sum(abbavec)
    baba = np.sum(babavec)
    if (abba + baba > 0):
        Dstat = (abba - baba) / (abba + baba)
    else:
        Dstat = float('nan')
    #fD
    checkfd1 = (p3_freq > p1_freq)
    abbafd1 = (1.0 - p2_freq)*p3_freq*p3_freq
    babafd1 = p2_freq*(1.0 - p3_freq)*p3_freq
    checkfd2 = (p3_freq < p1_freq)
    abbafd2 = (1.0 - p2_freq)*p1_freq*p1_freq
    babafd2 = p2_freq*(1.0 - p1_freq)*p1_freq
    abbafd = checkfd1 * abbafd1 + checkfd2 * abbafd2
    babafd = checkfd1 * babafd1 + checkfd2 * babafd2
    abbafd = np.sum(abbafd)
    babafd = np.sum(babafd)
    if (abbafd + babafd > 0):
        fD = (abba - baba) / (abbafd - babafd)
    else:
        fD = float('nan')
    #Heterozygosity
    hetvec = 2 * p3_freq * (1.0 - p3_freq)
    Het = np.sum(hetvec) / len_genome
    #Divergence
    #div = divergence ratio between pairs of individual sequences (admixed/source), 
    #and average across all pairs of individuals
    divratio = []
    for archi in range(0, p1_hap.shape[0]): #iterate over 0-99 haps; 100 total
        #divratio_check = []
        #print(archi)       
        divarchintro = vSumFunc(p3_hap, archi)
        divarchintro = divarchintro.astype("float")
        divarchnonintro = vSumFunc(p2_hap, archi)        
        divarchnonintro = divarchnonintro.astype("float") #took the inversion here so that the multiplying below is really actually dividing; probably not necessary       
        #pairwise combos of divarchintro/divarchnonintro -> 100*100 comparisons
        for comb in itertools.product(divarchintro,divarchnonintro): #pairwise combos of divarchintro/divarchnonintro
            #print(comb)
            if comb[1] != 0:
                #print(comb)
                divratio.append(comb[0]/comb[1])
                #divratio_check.append(comb[0]/comb[1])
        #print(float(sum(divratio_check)) / float(len(divratio_check)))        
        #[x if (x != float('inf')) for x in divratio]
    divratioavg = float(sum(divratio)) / float(len(divratio)) #len(divratio) = 100* (100*100)
    #print(divratioavg)
    #Q Dist
    Arc100 = (p1_freq == 1)
    NonAdm1 = (p2_freq < 0.01) 
    NonAdm10 = (p2_freq < 0.1)
    Arc100NonAdm1 = (Arc100 & NonAdm1)
    Arc100NonAdm10 = (Arc100 & NonAdm10)
    Freqs_Arc100NonAdm10 = p3_freq[np.where(Arc100NonAdm10 == True)]
    Freqs_Arc100NonAdm1 = p3_freq[np.where(Arc100NonAdm1 == True)]
    if Freqs_Arc100NonAdm10.size > 0:
        Q_10_100_q95 = np.percentile(Freqs_Arc100NonAdm10,95)
        Q_10_100_q90 = np.percentile(Freqs_Arc100NonAdm10,90)
        Q_10_100_max = np.max(Freqs_Arc100NonAdm10)
    else:
        Q_10_100_q95 = float('nan')
        Q_10_100_q90 = float('nan')
        Q_10_100_max = float('nan')
    if Freqs_Arc100NonAdm1.size > 0:
        Q_1_100_q95 = np.percentile(Freqs_Arc100NonAdm1,95)
        Q_1_100_q90 = np.percentile(Freqs_Arc100NonAdm1,90)
        Q_1_100_max = np.max(Freqs_Arc100NonAdm1)
    else:
        Q_1_100_q95 = float('nan')
        Q_1_100_q90 = float('nan')
        Q_1_100_max = float('nan')
    #Uniquely Shared Alleles
    U_10_0_100 = ( Arc100NonAdm10 & (p3_freq > 0) )
    U_10_20_100 = ( Arc100NonAdm10 & (p3_freq > 0.2) )
    U_10_50_100 = ( Arc100NonAdm10 & (p3_freq > 0.5) )
    U_10_80_100 = ( Arc100NonAdm10 & (p3_freq > 0.8) )		
    U_1_0_100 = ( Arc100NonAdm1 & (p3_freq > 0) )
    U_1_20_100 = ( Arc100NonAdm1 & (p3_freq > 0.2) )
    U_1_50_100 = ( Arc100NonAdm1 & (p3_freq > 0.5) )
    U_1_80_100 = ( Arc100NonAdm1 & (p3_freq > 0.8) )		
    U_10_0_100 = np.sum(U_10_0_100)
    U_10_20_100 = np.sum(U_10_20_100)
    U_10_50_100 = np.sum(U_10_50_100)
    U_10_80_100 = np.sum(U_10_80_100)        
    U_1_0_100 = np.sum(U_1_0_100)
    U_1_20_100 = np.sum(U_1_20_100)
    U_1_50_100 = np.sum(U_1_50_100)
    U_1_80_100 = np.sum(U_1_80_100)
    return Dstat, fD, Het, divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100, U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100



#######################################################
#additional stats at population level
def numseg (pop_hap):
    pop_freq = calc_freq (pop_hap)    
    pop_S = np.sum((pop_freq>0)&(pop_freq<1))    
    S_freq = pop_freq[(pop_freq>0)&(pop_freq<1)]    
    return pop_S,S_freq 

def thetaW (numS, popsize):
    x = list(range(1, popsize*2))
    denom = [1/float(i) for i in x]    
    thetaW = numS/sum(denom)
    return thetaW

def thetapi (S_freq, popsize):
    pi = sum(2*S_freq*(1.-S_freq))   
    n=popsize *2    
    thetapi = pi*n/(n-1)
    return thetapi

def thetaH (pophap, popsize):
    hap = np.transpose(pophap)
    n= popsize*2
    Si = [sum(i) for i in hap]
    Si = np.asarray([sum(i) for i in hap])
    num = [x**2*sum(Si == x) for x in list(range(1,n))]
    thetaH = sum(num)*2/(n*(n-1))
    return thetaH

    
def hap_homo (pophap,popsize):
    n=popsize*2   
    unique,count = np.unique(pophap, return_counts=True,axis=0)
    if len(unique) > 1:        
        freq = count/n
        freq[::-1].sort()
        p1 = freq[0]
        p2 = freq[1]
        H1 = sum(freq**2)
        H2 = H1-p1**2
        H12 = H1+2*p1*p2
        H2H1 = H2/H1    
    else:H1,H2,H12,H2H1=1,0,1,0
    return H1,H2,H12,H2H1

def stats_segthetasH12 (p1_hap,p2_hap,p3_hap, p1size,p2size,p3size):
    p1_S, p1_Sfreq= numseg(p1_hap)
    p2_S, p2_Sfreq= numseg(p2_hap)
    p3_S, p3_Sfreq= numseg(p3_hap)
    thetaW1, thetaW2, thetaW3 = thetaW (p1_S, p1size),thetaW (p2_S, p2size),thetaW (p3_S, p3size)
    thetapi1, thetapi2, thetapi3 = thetapi(p1_Sfreq,p1size),thetapi(p2_Sfreq,p2size),thetapi(p3_Sfreq,p3size)
    thetaH1, thetaH2, thetaH3 = thetaH (p1_hap, p1size),thetaH (p2_hap, p2size),thetaH (p3_hap, p3size)
    p1_H1, p1_H2,p1_H12, p1_H2H1 =  hap_homo (p1_hap, p1size)     
    p2_H1, p2_H2,p2_H12, p2_H2H1 =  hap_homo (p2_hap, p2size)
    p3_H1, p3_H2,p3_H12, p3_H2H1 =  hap_homo (p3_hap, p3size)         
    return p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1
    #24 total

def r2(pAB,paB,pAb,pab):
    pA = (pAB+pAb)
    pB = (pAB+paB)
    pa = (paB+pab)
    pb = (pAb+pab)
    D = pAB-pA*pB
    Dab = pab-pa*pb
    if (pA*(1-pA)*pB*(1-pB) != 0 ):
        r2 = (D**2)/(pA*(1-pA)*pB*(1-pB))
    else:
        r2= (Dab**2)/(pa*(1-pa)*pb*(1-pb))
    return r2

def Zns (pophap,popgeno):
    #popgeno = geno1
    pop_freq = calc_freq (pophap)    
    geno = np.transpose(popgeno)
    #S_hap = np.asarray([hap[(pop_freq>0)&(pop_freq<1)] for hap in pophap])
    S_geno = np.asarray([g[(pop_freq>0)&(pop_freq<1)] for g in geno])    
    #[g[41] for g in S_geno]    
    S = len(S_geno[0])
    if S>1:
        ls = range(1,S+1)
        pairs = list(itertools.combinations(ls,2))    
        geno_list = np.asarray([[g[n].replace('|', '') for n in range(0,len(S_geno[0]))] for g in S_geno])
        genos_combos = [[list(zip(g[x-1],g[y-1])) for g in geno_list] for x,y in pairs]    
        geno = [[''.join(i) for n in g for i in n] for g in genos_combos]
        AB = [sum([y == '11' for y in x]) for x in geno]
        aB = [sum([y == '01' for y in x]) for x in geno]
        Ab = [sum([y == '10' for y in x]) for x in geno]
        ab = [sum([y == '00' for y in x]) for x in geno]
        ps = [[w/float((w+x+y+z)),x/float((w+x+y+z)),y/float((w+x+y+z)),z/float((w+x+y+z))] for w,x,y,z in list(zip(AB,aB,Ab,ab))]
        #for n in range(len(ps)):
            #w,x,y,z=ps[n]
            #print(r2(w,x,y,z))
        r2s = [r2(w,x,y,z) for w,x,y,z in ps]
        Zns = sum(r2s)*2./(S*(S-1))
    else: 
        Zns = "Nan"    
    return Zns


def find_mut_site (segfile):
    segs = open(segfile)
    starts = []
    ends = []
    total = 0
    for line_counter,line in enumerate (segs):
        if line[0:4]=="exon":      
            fields = line.split()
            if (int(fields[1]) >= 2200000) & (int(fields[2]) <=2800000):
                starts.append(fields[1])
                ends.append(fields[2])
                total+=1        
    #mid = int(total/2)
    any_exon = random.choice(range(0,total))
    window_start = int(starts[any_exon])
    window_end = int(ends[any_exon])
    segs.close()
    return window_start,window_end


def find_mut_site_1 (segfile):
    segs = open(segfile)
    starts = []
    ends = []
    total = 0
    for line_counter,line in enumerate (segs):
        if line[0:4]=="exon":      
            fields = line.split()
            if (int(fields[1]) >= 2200000) & (int(fields[2]) <=2800000):
                starts.append(fields[1])
                ends.append(fields[2])
                total+=1        
    mid = int(total/2)
    window_start = int(starts[mid-1])
    window_end = int(ends[mid-1])
    segs.close()
    return window_start,window_end

#window_start,window_end=find_mut_site ("segments/sim_seq_info_bnc2.txt")

#temp_par = "sim_cluster_additive_adaptive.slim"
#region_name = "bnc2"
#new_par = "test_bnc2.slim"
    #update_par_file(region_name, temp_par, new_par)
def update_par_file(region_name, temp_par, new_par,insert,dominance,adm_amount,adm_time,sel_time): #dominance is a string = "additive" or "recessive" or "partial"
    if (dominance == "recessive")|(dominance == "h"):
        readLine = 18
        outtempLine = 62
        insertLine = 95
        readtempLine = 108
        #insertLine2 = 90
        treeLine = 173
        admLine = 133
        adm_timeLine = 130
        endadm_timeLine = 136
        sel_timeLine1 = 66
        sel_timeLine2 = 75
        sel_timeLine3 = 153
    else:
        readLine = 18
        outtempLine = 70
        insertLine = 103
        readtempLine = 116
        #insertLine2 = 98
        treeLine = 181
        admLine = 141
        adm_timeLine = 138
        endadm_timeLine = 144
        sel_timeLine1 = 74
        sel_timeLine2 = 83
        sel_timeLine3 = 161
    oldfile = open(temp_par)
    newfile = open(new_par,'w')
    line_counter=0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()
        if line_counter==readLine:
            #print(fields)
            fields[2] = 'readFile("/u/scratch/x/xinjunzh/slim_nonAIsweep/1ksegments/sim_seq_info_'+str(region_name)+'.txt");'
        elif line_counter == outtempLine:
            fields[0] = 'sim.outputFull("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_temp/'+str(region_name)+'_" '
            fields[1] = '+"misspec-back2afr_"+'
        elif line_counter == insertLine:
            fields[1] = str(insert)+");"
        elif line_counter == readtempLine:
            fields[0] = 'sim.readFromPopulationFile("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_temp/'+str(region_name)+'_" '
            fields[1] = '+"misspec-original_"+'          
        elif line_counter == treeLine:
            fields[0] = 'sim.treeSeqOutput("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_tree/'+str(region_name)+'_" '
            fields[1] = '+"misspec-original_"+'          
        elif line_counter == adm_timeLine:
            fields[0] = str(adm_time) 
        elif line_counter == admLine:
            fields[2] = str(adm_amount)+ "," 
            fields[3] = str(1-adm_amount)+"));"
        elif line_counter == endadm_timeLine:
            fields[0] = str(adm_time+1)
        elif line_counter == sel_timeLine1:
            fields[1] = str(adm_time+sel_time)
            fields[0] = str(adm_time)+":"
        elif line_counter == sel_timeLine2:
            fields[0] = str(adm_time+sel_time+1)#+":"
        elif line_counter == sel_timeLine3:
            fields[0] = str(adm_time+sel_time+1)                                    
        new_line=str()
        for item in fields:
            new_line = new_line+item+" "
        newfile.write(new_line+'\n')   
    newfile.close()
    oldfile.close()

def update_par_back2afr_file(region_name, temp_par, new_par,insert,dominance,adm_amount,adm_time,sel_time,back2afr_amount): #dominance is a string = "additive" or "recessive" or "partial"
    if (dominance == "")|(dominance == "h"):
        readLine = 18
        outtempLine = 62
        insertLine = 95
        readtempLine = 108
        #insertLine2 = 90
        treeLine = 173
        admLine = 133
        adm_timeLine = 130
        endadm_timeLine = 136
        sel_timeLine1 = 66
        sel_timeLine2 = 75
        sel_timeLine3 = 153
    else:
        readLine = 18
        outtempLine = 70
        insertLine = 103
        readtempLine = 116
        #insertLine2 = 98
        treeLine = 181#181
        admLine = 141
        #admLine2 = 150 #add second pulse at fixed time
        adm_timeLine = 138
        #adm_timeLine2 = 149
        endadm_timeLine = 144
        backmigration_Line = 166
        sel_timeLine1 = 74
        sel_timeLine2 = 83
        sel_timeLine3 = 161
    oldfile = open(temp_par)
    newfile = open(new_par,'w')
    line_counter=0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()
        if line_counter==readLine:
            #print(fields)
            fields[2] = 'readFile("/u/scratch/x/xinjunzh/slim_nonAIsweep/1ksegments/sim_seq_info_'+str(region_name)+'.txt");'
        elif line_counter == outtempLine:
            fields[0] = 'sim.outputFull("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_temp/'+str(region_name)+'_" '
            fields[1] = '+"misspec-back2afr_"+'
        elif line_counter == insertLine:
            fields[1] = str(insert)+");"
        elif line_counter == readtempLine:
            fields[0] = 'sim.readFromPopulationFile("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_temp/'+str(region_name)+'_" '
            fields[1] = '+"misspec-back2afr_"+'          
        elif line_counter == treeLine:
            fields[0] = 'sim.treeSeqOutput("/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/'+dominance+'/test_tree/'+str(region_name)+'_" '
            fields[1] = '+"misspec-back2afr_"+'          
        elif line_counter == adm_timeLine:
            fields[0] = str(adm_time) 
        elif line_counter == admLine:
            fields[2] = str(adm_amount)+ "," 
            fields[3] = str(1-adm_amount)+"));"
        elif line_counter == endadm_timeLine:
            fields[0] = str(adm_time+1)
#second pulse
        elif line_counter == backmigration_Line:
            fields[1] = "c("+str(back2afr_amount)+"));"   
        elif line_counter == sel_timeLine1:
            fields[1] = str(adm_time+sel_time)
            fields[0] = str(adm_time)+":"
        elif line_counter == sel_timeLine2:
            fields[0] = str(adm_time+sel_time+1)#+":"
        elif line_counter == sel_timeLine3:
            fields[0] = str(adm_time+sel_time+1)                                    
        new_line=str()
        for item in fields:
            new_line = new_line+item+" "
        newfile.write(new_line+'\n')   
    newfile.close()
    oldfile.close()


#######################################################

if __name__ == "__main__":
    os.chdir("/u/scratch/x/xinjunzh/slim_nonAIsweep/")    
    if whichh==1:
    	dominance="additive"
    elif whichh==2:
    	dominance="partial"
    elif whichh==3:
    	dominance="recessive"    
    #os.chdir("/Users/xinjunzhang/Dropbox/extract_s/")
    #thisfile="/Users/xinjunzhang/Desktop/AI-ML_Project/adaptiveIntrogressionML-master/scripts_simulation_treeseq/extract_s/segments/sim_seq_info_chr3region.txt"
    segs_path = "1ksegments/"
    allfiles = [f for f in os.listdir(segs_path) if f.endswith("txt")]
    genes = [f.partition("info_")[2].partition(".txt")[0] for f in allfiles] #len(genes) = 26
    thisgene = genes[whichgene-1]
    thisfile = allfiles[whichgene-1]
    DIR_stat = "/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+dominance+"/test_stat/"
    num_rep = 5    
    #slim_template_local = "/Users/xinjunzhang/Desktop/AI-ML_Project/adaptiveIntrogressionML-master/scripts_simulation_treeseq/extract_s/1kseg_exonrdist/misspec/misspec-back2afr_additive_adaptive.slim"
    #slim_thisgene_local = "/Users/xinjunzhang/Desktop/AI-ML_Project/adaptiveIntrogressionML-master/scripts_simulation_treeseq/extract_s/1kseg_exonrdist/misspec/misspec-back2afr_additive_test0.slim"    
    slim_template = "/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+"misspec-back2afr_"+dominance+"_adaptive.slim"
    slim_thisgene = "/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+dominance+"/"+thisgene+"-misspec-back2afr_"+dominance+"_adaptive_newtest.slim"    
    adm_range = [0.01,0.02,0.05,0.1]
    back2afr_range=[1e-5,5e-5,1e-4,5e-4,1e-3]
    #time-adm_range = range(0,100)
    rep=0
    while rep < num_rep:
        print(rep)
        try:
            back2afr_amount = back2afr_range[rep]
            adm_amount = random.choice(adm_range)
            window_start,window_end=find_mut_site ("1ksegments/"+thisfile)
            insert = int((int(window_end)+int(window_start))/2)   
            adm_time = random.choice(range(8697,8747))   
            sel_time = random.choice(range(0,61))   #up to 15000 years 
            update_par_back2afr_file(thisgene, slim_template, slim_thisgene,insert,dominance,adm_amount,adm_time,sel_time,back2afr_amount)
#1. run slim, get tree file
            slim_output = "/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+dominance+"/"+thisgene +"_misspec-back2afr_"+ "-slim_output_"+dominance+"_adaptive.txt"
            os.system('/u/home/x/xinjunzh/slim_build/slim %s > %s' %(slim_thisgene,slim_output))
            print("slim done")
            m2_s,m2_in_p4 = getSelectedAlleleSelcoeff(slim_output)
            m2_s = float(m2_s)

#2. from tree to vcf
            #filename = "/Users/xinjunzhang/Desktop/chr3region_misspec-back2afr_test_additive_adaptive"
            filename = "/u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+dominance+"/test_tree/"+thisgene +"_misspec-back2afr_"+ 'test_'+dominance+'_adaptive'
            scalingfactor = 10
            mu = 1.5e-8 * scalingfactor
            archaic_sample_time = 152.0

    #number of individuals from [afr, asn]
            sampsize=[108,99]#[100, 100] #Neanderthal is 1

    #get exon definitions
            #exonfilepath=thisfile
            exonfilepath = segs_path+thisfile #'sim_seq_info_EPAS1.txt'
            lines = open(exonfilepath, 'r').readlines()
            lines = [x for x in lines if x.startswith('exon')]
            lines = [x.strip('\n').split(' ') for x in lines]
            annot, start, end = zip(*lines)
            start = [int(x) for x in start]
            end = [int(x) for x in end]

            ts_sample = simplify_tree('{0}.trees'.format(filename), sampsize, archaic_sample_time)

            ns_vcf(ts_sample, '{0}_ns.vcf'.format(filename))


            noncoding_vcf(ts_sample, '{0}_noncoding.vcf'.format(filename), mu, start, end)
            syn_vcf(ts_sample, '{0}_syn.vcf'.format(filename), mu, start, end)

            cmd = '''
            cat {0}_ns.vcf | grep "##" > {0}_unsorted.vcf;
            echo '##INFO=<ID=TT,Number=A,Type=String,Description="Annotation">' >> {0}_unsorted.vcf;
            cat {0}_ns.vcf | grep "#CHROM" >> {0}_unsorted.vcf;
            cat {0}_ns.vcf | grep -v "#" | awk '$8="TT=NS"' | tr ' ' '\t' >> {0}_unsorted.vcf;
            cat {0}_syn.vcf | grep -v "#" | awk '$8="TT=SYN"' | tr ' ' '\t' >> {0}_unsorted.vcf;
            cat {0}_noncoding.vcf | grep -v "#" | awk '$8="TT=NC"' | tr ' ' '\t' >> {0}_unsorted.vcf;
            bcftools sort {0}_unsorted.vcf -o {0}_combined.vcf
            '''.format(filename)

            os.system(cmd)

            print("vcf done")

#from vcf to stats 
            vcfpath = filename + "_combined.vcf"
            os.rename(vcfpath, filename +"_rep"+str(rep)+ "_combined.vcf")
            vcfpath=filename +"_rep"+str(rep)+ "_combined.vcf"
            os.system("cp "+vcfpath+" /u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/"+dominance+"/test_VCF/")
            
            
            treepath = filename+".trees"
            len_genome = 4999999
            p1size,p2size,p3size = 2,108,99#100,100

            meanp1 = calc_p1ancestry(treepath,2,4108,8900-adm_time-1)
    
            geno1, geno2,geno3,len_genome,mut_types, pos_list = vcf2genos(vcfpath,p1size,p2size,p3size)
    
            hapmat1 = geno2hap(geno1)
            hapmat2 = geno2hap(geno2)
            hapmat3 = geno2hap(geno3)   
        
            p1_freq_all = calc_freq (hapmat1)
            p3_freq_all = calc_freq (hapmat3)
            locus = pos_list.index(insert)
            if m2_in_p4 !=0:
                m2freq_p1 = p1_freq_all[locus]
                m2freq_p3 = p3_freq_all[locus]
            else:  
                m2freq_p1 = p1_freq_all[locus]
                m2freq_p3 = 0
    
    #define exon windows, subset hapmat
            start_w = 0
    #end_w = 0
            windows = []
            start_pos = []
    #each = []
            i=0
            count=0
            while i <= 5000000:#4999999:
                if (count%50000==1)& (start_w == 0):
                    start_w +=1
                    start_pos.append(i)
                    i+=1
                    count+=1
                if (count%50000==0) & (start_w > 0):
                    end_pos = i
                    if start_w > 0:
                        windows.append([start_pos[0],end_pos])
                #print(start_pos[0],end_pos)
                    start_w = 0
                    i=start_pos[0]+9999
                    start_pos = []
                    count = 0
                    #i+=1
        
                else:
                    i+=1
                    count+=1

    #len(windows)

            focal_window = [w for w in windows if insert in range(w[0],w[1]+1)][0]
            outputname = thisgene+"_misspec-back2afr_"+dominance+"_adaptive_"
            with open(DIR_stat+thisgene+"/"+str(rep)+outputname+batch+".txt", "w") as stats_file:
                stats_name = ["rep","classifier","thisgene","distance","adm_time","adm_amount","sel_time","m2_s","m2freq_p1","m2freq_p3","start","end","num_alleles","meanp1","Dstat","fD","Het",
            	"divratioavg","Q_10_100_q95","Q_10_100_q90","Q_10_100_max","Q_1_100_q95","Q_1_100_q90","Q_1_100_max",
            	"U_10_0_100","U_10_20_100","U_10_50_100","U_10_80_100","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100",
            	"p1_S","p2_S","p3_S","thetaW1","thetaW2","thetaW3","thetapi1","thetapi2","thetapi3","thetaH1","thetaH2","thetaH3",
            	"p1_H1","p1_H2","p1_H12","p1_H2H1","p2_H1","p2_H2","p2_H12","p2_H2H1","p3_H1","p3_H2","p3_H12","p3_H2H1","Zns1","Zns2","Zns3"]

                stats_file.writelines(i+"\t" for i in stats_name)
                stats_file.write("\n")
        
                for w in range(len(windows)):
                    thesepos = [pos for pos in pos_list if (pos <= windows[w][1] and pos >= windows[w][0])]
            
                    start = pos_list.index(thesepos[0])
                    end = pos_list.index(thesepos[len(thesepos)-1])
                    len_seg = windows[w][1]-windows[w][0]
            #stats_file.writelines(str(s)+"\t" for s in windows[w])
            

                    p1_hap,p2_hap,p3_hap = hapmat1[:,start:end],hapmat2[:,start:end],hapmat3[:,start:end]
            
                    Dstat,fD,Het,divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,U_10_80_100,U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100 = calc_stats(p1_hap,p2_hap,p3_hap,len_seg)

                    p1_S, p2_S, p3_S, thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1 = stats_segthetasH12 (p1_hap,p2_hap,p3_hap, p1size,p2size,p3size)

            
                    geno_p1,geno_p2,geno_p3 = geno1[start:end],geno2[start:end],geno3[start:end]
            
                    Zns1 = Zns (p1_hap,geno_p1)
                    Zns2 = Zns (p2_hap,geno_p2)
                    Zns3 = Zns (p3_hap,geno_p3)
                    dist = abs(((windows[w][0]+windows[w][1])/2)-insert)
                
                    classifier = 0
                    if insert in range(windows[w][0],windows[w][1]):
                        classifier = 1
                    elif (((windows[w][0]+windows[w][1])/2) >= insert-250000) & (windows[w][0] < focal_window[0]):
                        classifier = 2
                    elif (windows[w][0] > focal_window[0]) & (((windows[w][0]+windows[w][1])/2) <= insert+250000):
                        classifier = 2                
            
                    stats_list = [rep,classifier,thisgene,dist,adm_time,adm_amount,sel_time,m2_s,m2freq_p1,m2freq_p3,windows[w][0],windows[w][1],len(thesepos),meanp1,Dstat,fD,Het,divratioavg,Q_10_100_q95,Q_10_100_q90,Q_10_100_max,
                              Q_1_100_q95,Q_1_100_q90,Q_1_100_max,U_10_0_100,U_10_20_100,U_10_50_100,
                              U_10_80_100,U_1_0_100,U_1_20_100,U_1_50_100,U_1_80_100,p1_S, p2_S, p3_S, 
                              thetaW1, thetaW2, thetaW3,thetapi1, thetapi2, thetapi3,thetaH1, thetaH2, thetaH3,
                              p1_H1, p1_H2,p1_H12, p1_H2H1,p2_H1, p2_H2,p2_H12, p2_H2H1,p3_H1, p3_H2,p3_H12, p3_H2H1,Zns1,Zns2,Zns3]

            
                    stats_file.writelines(str(i)+"\t" for i in stats_list)
                    stats_file.write("\n")
                    stats_file.flush()
            print("stats done")
            rep+=1
        	
        except:
            pass
    print("END OF SIMULATION")
    #p1_hap,p2_hap,p3_hap = hapmat1,hapmat2,hapmat3













