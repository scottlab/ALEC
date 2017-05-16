#!/usr/bin/env python
#title               :ALEC.py
#description         :Correct long read sequencing(Pacbio or Nanopore)
#author              :Yao Yang
#date of 1st version :05/30/2015
#last revised date   :05/10/2017
#version             :1.0
#usage               :python ALEC.py
#notes               :
#python_version      :2.7.10
#================================================================================================================================================================================

'''Import the modules needed to run the script'''
import sys
import re
import time
import pysam
import argparse
import subprocess
import math
import numpy
import random
from itertools import groupby
from operator import itemgetter
from itertools import product
start_time = time.time()

def get_argument():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required = True,
                        help='Input file, currently SAM only, but will expand to BAM')
    parser.add_argument('-r', '--reference', required = True,
                        help='Reference file')
    parser.add_argument('-t', '--targetRegion', required = True,
                        help='Target regions interval')
    parser.add_argument('-del', '--deletion',  
                        help='Deletion error frequency baseline to trigger correction.', 
                        type = float, default = 0.15)
    parser.add_argument('-ins', '--insert', 
                        help='Insert error frequency baseline to trigger correction.', 
                        type = float, default = 0.15)
    parser.add_argument('-mis', '--mismatch', 
                        help='Mismatch error frequency baseline to trigger correction.', 
                        type = float, default = 0.05)
    parser.add_argument('-lf', '--lengthFilter', 
                        help='Reads shorter than lengthFilter*Reference_length will be excluede in this correction process', 
                        type = float, default = 0.0)
    parser.add_argument('-del_hp', '--del_homo_p', 
                        help='Deletion Homopolymer Penatly.', 
                        type = float, default = 0.05)
    parser.add_argument('-ins_hp', '--ins_homo_p', 
                        help='Insert Homopolymer Penatly.', 
                        type = float, default = 0.0)
    parser.add_argument('-ds', '--downsample_freq', 
                        help='downsampling_freq', 
                        type = float, default = 1.0)
    parser.add_argument('-x', '--platform',choices=['pacbio_ccs','pacbio_sub','nanopore'],
                        help='Platform and type of sequencing data')
    args = parser.parse_args()
    if args.platform == 'pacbio_ccs':
        [args.deletion, args.insert, args.mismatch, args.del_homo_p]  = [0.12,0.1,0.1,0.12]
    elif args.platform == 'pacbio_sub':
        [args.deletion, args.insert, args.mismatch, args.del_homo_p]  = [0.18,0.15,0.1,0.09]
    elif args.platform == 'nanopore':
        [args.deletion, args.insert, args.mismatch, args.del_homo_p]  = [0.30,0.2,0.2,0.06] 

def parse_fasta(reference_file):
    '''getting DNA sequence from reference fasta file'''
    global ref_seq, len_ref, ref_context_info,region_chr, region_start, region_end
    counter = 1
    region_chr = args.targetRegion[:args.targetRegion.find(":")]
    region_start = int(args.targetRegion[args.targetRegion.find(":")+1:args.targetRegion.find("-")])
    region_end = int(args.targetRegion[args.targetRegion.find("-")+1:])
    ref_context_info = []
    referencefile = pysam.Fastafile(args.reference)         
    ref_seq = referencefile.fetch(region = args.targetRegion)
    for i in range(len(ref_seq)-1):
        if ref_seq[i] == ref_seq[i+1]:
            counter += 1
            if i == len(ref_seq)-2:
                for j in range(counter):
                    ref_context_info.append([ref_seq[i], counter, counter-j])
        else:
            for j in range(counter):
                ref_context_info.append([ref_seq[i], counter,counter-j])            
            counter = 1
            if i == len(ref_seq)-2:
                ref_context_info.append([ref_seq[-1],1,1])
    len_ref = region_end-region_start+1
    return [ref_seq,ref_context_info]

def RevC(seq):
    '''getting DNA reverse complement'''
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N','D':'D'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def parse_sam(read):
    '''extracting single read information from each line in sam file'''
    match, read_seq  = (str() for i in range(2))
    base_length,operator,insert_bases = ([] for i in range(3))
    insert = False
    read_seq = str()
    cigar_tuple = read.cigartuples
    seq_trimmed = read.query_alignment_sequence
    for i in range(len(cigar_tuple)):
        if cigar_tuple[i][0] == 0:
            read_seq += seq_trimmed[:cigar_tuple[i][1]]
            seq_trimmed = seq_trimmed[cigar_tuple[i][1]:]
            if not insert:
                match += '1'*(cigar_tuple[i][1])
                insert_bases += ['O']*(cigar_tuple[i][1]) 
            else:
                match += '1'*(cigar_tuple[i][1]-1)
                insert_bases += ['O']*(cigar_tuple[i][1] - 1)
                insert = False
        elif cigar_tuple[i][0] == 2:
            read_seq += 'D'*cigar_tuple[i][1]
            match += '0'*cigar_tuple[i][1]
            insert_bases += ['']*(cigar_tuple[i][1])
        elif cigar_tuple[i][0] == 1:
            match += '2'
            insert_bases.append(seq_trimmed[:cigar_tuple[i][1]])
            seq_trimmed = seq_trimmed[cigar_tuple[i][1]:]
            insert = True 
    if read.reference_start+1 > region_start:
        match = '.' * (read.reference_start + 1 -region_start) + match
        insert_bases = ['']*(read.reference_start + 1 -region_start)+ insert_bases
        read_seq = '.'* (read.reference_start + 1 -region_start) +read_seq 
    else:
        match = match[(region_start-read.reference_start-1):] 
        insert_bases = insert_bases[(region_start-read.reference_start-1):]
        read_seq = read_seq[(region_start-read.reference_start-1):]
    match = match[:len_ref]
    insert_bases = insert_bases[:len_ref]
    match += '.'*(len_ref - len(match))
    insert_bases += [''] * (len_ref - len(insert_bases))
    read_seq = read_seq[:len_ref]
    read_seq +=  '.'*(len_ref - len(read_seq))
    return [match,read_seq,insert_bases]

def extract_feature_matrix(reads_in_sam):
    '''extract match_maatrix, base_matrix and insert_matrix from sam file'''
    global base_matrix, match_matrix,insert_bases_matrix, read_ct 
    match_matrix,base_matrix,insert_bases_matrix,reads_ID,flag  = ([] for i in range(5))
    for line in reads_in_sam:
        if (line.reference_length > args.lengthFilter * len_ref
            and samfile.getrname(line.reference_id) == region_chr
            and line.reference_start < region_end
            and line.reference_end > region_start):
            parsed_read = parse_sam(line)
            match_matrix.append(parsed_read[0])
            base_matrix.append(parsed_read[1])
            insert_bases_matrix.append(parsed_read[2])
            reads_ID.append(line.query_name)
            flag.append(line.flag)
    match_matrix = map(list,zip(*match_matrix))
    base_matrix = map(list,zip(*base_matrix))
    insert_bases_matrix = map(list,zip(*insert_bases_matrix))
    read_ct = len(flag)
    return [match_matrix, base_matrix, insert_bases_matrix,reads_ID,flag]

def get_allele_freq_lst(base_matrix):
    global ct_base_lst, consensus_base, depth,ct_mismatch,ct_mismatch_lst
    ct_base_lst,allele_freq_lst,depth,ct_mismatch_lst,consensus_base = ([] for i in range(5))
    for i in range(len_ref):
        ct_mismatch=0
        A_ct = base_matrix[i].count('A')
        T_ct = base_matrix[i].count('T')
        C_ct = base_matrix[i].count('C')
        G_ct = base_matrix[i].count('G')
        D_ct = base_matrix[i].count('D')
        total_ct = float(A_ct + T_ct + C_ct + G_ct + D_ct)
        depth.append(total_ct)
        allele_freq = [A_ct/total_ct,T_ct/total_ct,C_ct/total_ct,G_ct/total_ct,D_ct/total_ct]
        ct_base_lst.append([A_ct,T_ct,C_ct,G_ct])
        for j in range(4):
            if allele_freq[j] < 0.1:
                ct_mismatch += [A_ct,T_ct,C_ct,G_ct][j]
        allele_freq_lst.append(allele_freq)
        ct_mismatch_lst.append(ct_mismatch)
        max_af = max(allele_freq)
        idx = allele_freq.index(max_af)
        consensus_base.append('ATCGD'[idx])
    return [allele_freq_lst,consensus_base,depth]

def get_indel_rate(match_matrix):
    global del_rate, ins_rate,ct_del_lst, ct_ins_lst,  ct_dep_lst
    ct_del_lst, ct_ins_lst,  ct_dep_lst,del_rate, ins_rate = ([] for i in range(5))
    for i in range(len_ref):
        match_count = match_matrix[i].count('1')
        del_count = match_matrix[i].count('0')
        ins_count = match_matrix[i].count('2')
        base_depth = match_count + del_count + ins_count
        del_rate.append(del_count/float(base_depth))
        ins_rate.append(ins_count/float(base_depth))
        ct_del_lst.append(del_count)
        ct_ins_lst.append(ins_count)
        ct_dep_lst.append(base_depth)
    return [del_rate,ins_rate]

def get_consensus_insert_bases(insert_bases_matrix):
    consensus_insert,sninsert_freq,A_insert_freq,T_insert_freq,C_insert_freq,G_insert_freq =( [] for i in range(6))
    sninsert_freq = []
    for i in range(len_ref):
        inserts_at_pos = list(set(insert_bases_matrix[i]))
        if '' in inserts_at_pos:
            inserts_at_pos.remove('')
        ct_insert_bases = []
        for element in inserts_at_pos:
            ct_insert_bases.append(insert_bases_matrix[i].count(element))
        if ct_insert_bases == []:
            consensus_insert.append('')
        elif inserts_at_pos[ct_insert_bases.index(max(ct_insert_bases))] != 'O':
            consensus_insert.append('')
        else:
            consensus_insert.append(inserts_at_pos[ct_insert_bases.index(max(ct_insert_bases))])
        A_insert_freq.append(insert_bases_matrix[i].count('A')/depth[i])
        T_insert_freq.append(insert_bases_matrix[i].count('T')/depth[i])
        C_insert_freq.append(insert_bases_matrix[i].count('C')/depth[i])
        G_insert_freq.append(insert_bases_matrix[i].count('G')/depth[i])
        
    return [consensus_insert,[A_insert_freq,T_insert_freq,C_insert_freq,G_insert_freq]]

def max_count_of_base(DNA_string):
    max_ct = 0
    nucleotide_set = []
    for i in DNA_string:
        if i not in nucleotide_set:
            nucleotide_set.append(i)
    for i in list(nucleotide_set):
        if DNA_string.count(i) > max_ct:
            max_ct = DNA_string.count(i)
            max_base = i
    return [max_ct,max_base]

def find_long_del():
    global true_large_del
    connected_deletion_matrix = [['0' for i in range(read_ct)] for j in range(len_ref-1)] 
    connected_deletion_freq =[0.0 for i in range(len_ref-1)]
    large_del_pos = []
    large_del_freq = []
    large_del = []
    true_large_del = []
    for i in range(len_ref-1):
        for j in range(read_ct):
            if match_matrix[i][j] == '0' and match_matrix[i+1][j] == '0':
                connected_deletion_matrix[i][j] = '1'
            if match_matrix[i][j] == '.' or match_matrix[i+1][j] == '.':
                connected_deletion_matrix[i][j] = '.'
    for i in range(len_ref-1):
        connected_deletion_freq[i] = float(connected_deletion_matrix[i].count('1'))/(connected_deletion_matrix[i].count('1')+connected_deletion_matrix[i].count('0'))
    for i in range(len(connected_deletion_freq)):
        if connected_deletion_freq[i] >= args.deletion:
            large_del_pos.append(i)
            large_del_freq.append(connected_deletion_freq[i]) 
    for k, g in groupby(enumerate(large_del_pos), lambda (i, x): i-x): 
        large_del.append(map(itemgetter(1), g))
    for pos in large_del:
        pos+=[pos[-1]+1]
    for i in range(len(large_del)):
        tem_del_freq = [0.0 for k in range(len(large_del[i]))]
        for j in range(len(large_del[i])):
            tem_del_freq[j] = connected_deletion_freq[large_del[i][j]]
        if len(tem_del_freq)>2 and numpy.std(tem_del_freq)<=0.1:
            true_large_del+=large_del[i]


def adjust_deletion_pos(i,j):
    '''correct deletion misplaced by aliner'''
    if (i>0 and base_matrix[i-1][j] not in ['.','R','D']
        and allele_freq_lst[i-1]['ATCG'.index(base_matrix[i-1][j])] <= args.mismatch
        and allele_freq_lst[i]['ATCG'.index(base_matrix[i-1][j])] >= args.mismatch):
        base_matrix[i][j]=base_matrix[i-1][j]
        match_matrix[i][j] = '1'
        base_matrix[i-1][j] = 'D'
        match_matrix[i-1][j] = '0'
    elif (i<len(base_matrix)-1
          and base_matrix[i+1][j] not in ['.','R','D']
          and allele_freq_lst[i+1]['ATCG'.index(base_matrix[i+1][j])] <= args.mismatch
          and allele_freq_lst[i]['ATCG'.index(base_matrix[i+1][j])] >= args.mismatch):
        base_matrix[i][j]=base_matrix[i-1][j]
        match_matrix[i][j] = '1'
        base_matrix[i+1][j] = 'D'
        match_matrix[i+1][j] = '0'

def revise_seq_context_info():
    for i in range(len_ref):
        if consensus_base[i] not in [ref_seq[i],'D']:
            if consensus_base[i] == ref_context_info[i+1][0]:
                ref_context_info[i][0] = ref_context_info[i+1][0]
                ref_context_info[i][1] += ref_context_info[i+1][1]
                ref_context_info[i][2] += ref_context_info[i+1][1]
            if consensus_base[i] == ref_context_info[i-1][0]:
                ref_context_info[i][0] = ref_context_info[i-1][0]
                ref_context_info[i][1] += ref_context_info[i-1][1]
                ref_context_info[i][2] += ref_context_info[i-1][1]

def replace_confusing_deletion(i,j):
    '''correct deletion misplaced by aliner'''
    del_len = 1
    while i+del_len< len_ref and match_matrix[i+del_len][j] == '0':
        del_len+=1
    if (i+2*del_len<len(base_matrix)
        and i+del_len not in true_large_del):
        for k in range(del_len):
            if base_matrix[i+k+del_len][j]!=consensus_base[i+k]:
                break
            else:
                base_matrix[i+k][j] = consensus_base[i+k]
                match_matrix[i+k][j] = '1'  
                base_matrix[i+k+del_len][j] = 'D'
                match_matrix[i+k+del_len][j] = '0'
 
def correct_del_in_new_2bp_homopolymer(i,j):
    '''correct deletion caused by homopolymer caseed by polymeorphism'''
    if (i>0 and base_matrix[i-1][j] not in ['.','R','D']
        and allele_freq_lst[i-1]['ATCG'.index(base_matrix[i-1][j])] >= args.mismatch
        and allele_freq_lst[i]['ATCG'.index(base_matrix[i-1][j])] >= args.mismatch):
        base_matrix[i][j]=base_matrix[i-1][j]
        match_matrix[i][j] = '1'
    elif (i<len(base_matrix)-1
          and base_matrix[i+1][j] not in ['.','R','D']
          and allele_freq_lst[i+1]['ATCG'.index(base_matrix[i+1][j])] >= args.mismatch
          and allele_freq_lst[i]['ATCG'.index(base_matrix[i+1][j])] >= args.mismatch):
        base_matrix[i][j]=base_matrix[i+1][j]
        match_matrix[i][j] = '1'

def correct_deletion(i):
    '''correct deletion error'''
    if i in true_large_del:
        pass
    elif del_rate[i] <=  args.deletion + args.del_homo_p*(ref_context_info[i][2]-1):
        if consensus_base[i] == 'D' and args.deletion + args.del_homo_p*(ref_context_info[i][2]-1) >=1:
            consensus_base[i] = ref_seq[i]
        for j in range(read_ct):
            if match_matrix[i][j] == '0':        
                adjust_deletion_pos(i,j)
            if match_matrix[i][j] == '0':
                replace_confusing_deletion(i,j)
            if match_matrix[i][j] == '0':
                base_matrix[i][j] = consensus_base[i]
                match_matrix[i][j] = '1'
    elif (i>0 
          and del_rate[i] >=  args.deletion + args.del_homo_p*(ref_context_info[i][2]-1)
          and del_rate[i] <= args.deletion + args.del_homo_p*ref_context_info[i][2]):
        for j in range(read_ct):
            if match_matrix[i][j] == '0':
                correct_del_in_new_2bp_homopolymer(i,j)
    elif del_rate[i] >= 1 - args.insert:
        for j in range(read_ct):
            base_matrix[i][j] = 'D'    
            match_matrix[i][j] = '0'

def correct_single_insert(i):
    for ins in ['A','T','C','G']:
        if snv_ins_freq['ATCG'.index(ins)][i]<= args.insert +args.ins_homo_p*(ref_context_info[i][2]-1):
            for j in range(read_ct):
                if insert_bases_matrix[i][j] == ins:
                    if (i>0 
                        and base_matrix[i-1][j] 
                        and allele_freq_lst[i-1]['ATCGD'.index(base_matrix[i-1][j])] <= args.mismatch 
                        and  insert_bases_matrix[i][j] == consensus_base[i-1]):
                        base_matrix[i-1][j]=consensus_base[i-1]
                        match_matrix[i-1][j] = '1'
                    elif (i<len(base_matrix)-1 
                          and base_matrix[i+1][j] 
                          and allele_freq_lst[i+1]['ATCGD'.index(base_matrix[i+1][j])] <= args.mismatch 
                          and  insert_bases_matrix[i][j] == consensus_base[i+1]):
                        base_matrix[i+1][j]=consensus_base[i+1]
                        match_matrix[i+1][j] = '1'
                    insert_bases_matrix[i][j] = 'O'
                    match_matrix[i][j] = '1'
   
def correct_large_insert(i):
    if ins_rate[i] <= args.insert +args.ins_homo_p*(ref_context_info[i][2]-1):
        for j in range(read_ct):
            if insert_bases_matrix[i][j] not in ['','O']:
                if (i>0 
                   and base_matrix[i-1][j] not in ['.','R'] 
                   and allele_freq_lst[i-1]['ATCGD'.index(base_matrix[i-1][j])] <= args.mismatch 
                   and  insert_bases_matrix[i][j] == consensus_base[i-1]):
                    base_matrix[i-1][j]=consensus_base[i-1]
                    match_matrix[i-1][j] = '1'
                elif (i<len(base_matrix)-1 
                      and base_matrix[i+1][j] not in ['.','R'] 
                      and allele_freq_lst[i+1]['ATCGD'.index(base_matrix[i+1][j])] <= args.mismatch 
                      and  insert_bases_matrix[i][j] == consensus_base[i+1]):
                    base_matrix[i+1][j]=consensus_base[i+1]
                    match_matrix[i+1][j] = '1'
                elif len(insert_bases_matrix[i][j]) == 1:
                    if (base_matrix[i][j] not in ['.','R'] 
                       and allele_freq_lst[i]['ATCGD'.index(insert_bases_matrix[i][j])] >= args.mismatch):
                        base_matrix[i][j]= insert_bases_matrix[i][j]
                insert_bases_matrix[i][j] = 'O'
                match_matrix[i][j] = '1'
    elif ins_rate[i] >= 1 - args.deletion:
        for j in range(read_ct):
            if insert_bases_matrix[i][j] == 'O':
                insert_bases_matrix[i][j] =consensus_insert[i]
                match_matrix[i][j] = '2'
  

def correct_miscall(match_matrix,consensus_base,base_matrix,allele_freq_lst):
    for i,j in product(range(len_ref),range(read_ct)):
        if  (base_matrix[i][j] not in ['.','R','D'] 
             and (allele_freq_lst[i]['ATCGD'.index(base_matrix[i][j])] <= args.mismatch)):
            if (i>0 and base_matrix[i-1][j] not in ['.','R','D']
                and allele_freq_lst[i-1]['ATCG'.index(base_matrix[i][j])] >= args.mismatch):
                base_matrix[i-1][j]=base_matrix[i][j]
            elif (i<len(base_matrix)-1
                and base_matrix[i+1][j] not in ['.','R','D']
                and allele_freq_lst[i+1]['ATCG'.index(base_matrix[i][j])] >= args.mismatch):
                base_matrix[i+1][j]=base_matrix[i][j]
            base_matrix[i][j] = consensus_base[i]
            match_matrix[i][j] = '1'
    return [base_matrix,match_matrix]  

def correct_reads(reads_in_sam):
    global mis_rate, allele_freq_lst, snv_ins_freq 
    mis_rate = []
    correction,allele_freq_lst  = ([] for i in range(2))
    match_matrix, base_matrix, insert_bases_matrix, reads_ID, flag = extract_feature_matrix(reads_in_sam)
    allele_freq_lst,consensus_base,depth = get_allele_freq_lst(base_matrix)
    [consensus_insert,snv_ins_freq] = get_consensus_insert_bases(insert_bases_matrix)
    for i in range(len_ref):
        mis_rate.append(1-allele_freq_lst[i]["ATCG".index(ref_seq[i])]-allele_freq_lst[i][4])
    del_rate,ins_rate = get_indel_rate(match_matrix)
    find_long_del()
    revise_seq_context_info() 
    for i in range(len_ref):
        correct_large_insert(i)
        correct_deletion(i)
    del_rate,ins_rate = get_indel_rate(match_matrix)
    for i in range(len_ref):
        correct_single_insert(i)
    del_rate,ins_rate = get_indel_rate(match_matrix)
    allele_freq_lst,consensus_base,depth = get_allele_freq_lst(base_matrix)
    [base_matrix,match_matrix] = correct_miscall(match_matrix,consensus_base,base_matrix,allele_freq_lst)
    allele_freq_lst,consensus_base,depth = get_allele_freq_lst(base_matrix)
    match_matrix = map(list,zip(*match_matrix)) 
    base_matrix = map(list,zip(*base_matrix))
    insert_bases_matrix = map(list,zip(*insert_bases_matrix))
    for i in range(read_ct):
        corrected_seq = ''
        for j in range(len_ref):
            if match_matrix[i][j] == '1':
                corrected_seq += base_matrix[i][j]
            elif match_matrix[i][j] == '2':
                corrected_seq += (insert_bases_matrix[i][j]+base_matrix[i][j])
        corrected_seq = corrected_seq.replace('D','')
        corrected_seq = corrected_seq.replace('R','')
        if flag[i] == 0 or flag[i] == 2048:
            correction.append(corrected_seq)
        else:
            correction.append(RevC(corrected_seq))
    return [reads_ID,correction]

def output_corrected_reads(cor_reads, input_file):
    with open(input_file[:-3]+"corrected.fasta","w+") as output_file:
        for i in range(len(cor_reads[1])):
            output_file.write('>' + cor_reads[0][i] + '\n' + cor_reads[1][i] + '\n')

def output_error_rate(input_file):
    with open(input_file[:-3]+"error_rate_file_after","w+") as error_rate_file:
        for i in range(len(mis_rate)):
            error_rate_file.write('\t'.join(map(str,[region_chr, i + region_start, round(mis_rate[i],5),round(del_rate[i],5),round(ins_rate[i],5)]+map(str,ct_base_lst[i])+map(str,[ct_mismatch_lst[i] , ct_del_lst[i],ct_ins_lst[i],ct_dep_lst[i]])+ ref_context_info[i]))+'\n')

def load_input_file():
    global samfile, reads
    samfile = pysam.AlignmentFile(args.input,"r")
    if args.input[-3:] in ['bam','BAM']:
        reads = samfile.fetch(region=args.targetRegion)
    else:
        reads = samfile.fetch()

def main():
    get_argument()
    print args
    print "#################################################################################"
    print "Correcting reads in region" + args.targetRegion
    parse_fasta(args.reference)
    load_input_file()
    time_b = time.time()
    print ("ALEC took %s seconds to get reference info" % ( round((time_b - start_time),0)))
    time_c = time.time()
    print ("ALEC took %s seconds to get reads" % ( round((time_c - time_b),0)))
    cor_reads = correct_reads(reads)
    output_corrected_reads(cor_reads, args.input)
    print ("ALEC took %s seconds to correct" % (round((time.time() - time_c),0)))
    print ("ALEC took %s seconds to correct %s reads in region %s" % ( round((time.time() - start_time),2), len(cor_reads[1]), args.targetRegion))
    print "#################################################################################"
if __name__ == "__main__":
    main()

