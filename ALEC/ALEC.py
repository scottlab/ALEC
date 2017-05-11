#!/usr/bin/env python
#title               :ALEC.py
#description         :Correct Long Read Sequencing(Pacbio or Nanopore)
#author              :Yao Yang
#date of 1st version :05/30/2015
#last revised date   :05/01/2017
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
start_time = time.time()

def get_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file, currently SAM only, but will expand to BAM")
    parser.add_argument("-r", "--reference", help="Reference file")
    parser.add_argument("-t", "--targetRegion", help="Target regions interval")
    parser.add_argument("-del", "--deletion", help="Deletion error frequency baseline to trigger correction.", type = float, default = 0.15)
    parser.add_argument("-ins", "--insert", help="Insert error frequency baseline to trigger correction.", type = float, default = 0.15)
    parser.add_argument("-mis", "--mismatch", help="Mismatch error frequency baseline to trigger correction.", type = float, default = 0.05)
    parser.add_argument("-lf", "--lengthFilter", help="Reads shorter than lengthFilter*Reference_length will be excluede in this correction process", type = float, default = 0.0)
    parser.add_argument("-del_hp", "--del_homo_p", help="Deletion Homopolymer Penatly.", type = float, default = 0.05)
    parser.add_argument("-ins_hp", "--ins_homo_p", help="Insert Homopolymer Penatly.", type = float, default = 0.0)
    parser.add_argument("-ds", "--downsample_freq", help="downsampling_freq", type = float, default = 1.0)
    return parser


def parse_fasta(reference_file):
    '''getting DNA sequence from reference fasta file'''
    global ref_seq, len_ref, ref_context_info
    counter = 1
    ref_seq = str()
    ref_context_info = []
    with open (str(region_chr) + "_" + str(region_start) + "_" + str(region_end) + ".fasta","w+") as  temp_target_file:
        cmd= ['samtools', 'faidx', args.reference, args.targetRegion]
        subprocess.call(cmd, stdout= temp_target_file)
        temp_target_file.seek(0)
        for line in temp_target_file:
            if line[0] != '>':
                line = line.rstrip()
                ref_seq += str(line)
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
    len_ref = len(ref_seq)
    return [ref_seq,ref_context_info]

def RevC(seq):
    '''getting DNA reverse complement'''
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N','D':'D'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def parse_sam(read):
    '''extracting single read information from each line in sam file'''
    read = read.rstrip()
    items = read.split('\t')
    [flag,start,cigar,seq] = [items[1],int(items[3]), items[5], items[9]]
    match, read_seq  = (str() for i in range(2))
    base_length,operator,insert_bases = ([] for i in range(3))
    insert = False
    read_seq = str()
    while len(cigar)>0:
        for i in range(len(cigar)):
            if cigar[i].isalpha():
                pos = i
                break
        base_length.append(int(cigar[:pos]))
        operator.append(cigar[pos])
        cigar = cigar[pos+1:]
    if operator[0] in ['S','H']:
        if operator[0] == 'S':
            seq = seq[base_length[0]:]
        operator.pop(0)
        base_length.pop(0)
    if operator[-1] in ['S','H']:
        if operator[-1] == 'S':
            seq = seq[:len(seq) - base_length[-1]]
        operator.pop()
        base_length.pop()
    seq_clipping_trimmed = seq
    for i in range(len(base_length)):
        if operator[i] == 'M':
            read_seq += seq_clipping_trimmed[:base_length[i]]
            seq_clipping_trimmed = seq_clipping_trimmed[base_length[i]:]
            if not insert:
                match += '1'*(base_length[i])
                insert_bases += ['O']*(base_length[i]) 
            else:
                match += '1'*(base_length[i]-1)
                insert_bases += ['O']*(base_length[i] - 1)
                insert = False
        elif operator[i] == 'D':
            read_seq += 'D'*base_length[i]
            match += '0'*base_length[i]
            insert_bases += ['']*(base_length[i])
        elif operator[i] == 'I':
            match += '2'
            insert_bases.append(seq_clipping_trimmed[:base_length[i]])
            seq_clipping_trimmed = seq_clipping_trimmed[base_length[i]:]
            insert = True 
    if start > region_start:
        match = '.' * (start-region_start) + match
        insert_bases = ['']*(start-region_start)+ insert_bases
        read_seq = '.'* (start-region_start) +read_seq 
    else:
        match = match[(region_start-start):] 
        insert_bases = insert_bases[(region_start-start):]
        read_seq = read_seq[(region_start-start):]
    match = match[:len_ref]
    insert_bases = insert_bases[:len_ref]
    match += '.'*(len_ref - len(match))
    insert_bases += [''] * (len_ref - len(insert_bases))
    read_seq = read_seq[:len_ref]
    read_seq +=  '.'*(len_ref - len(read_seq))
    read_length_overlapped_with_target_region = len([c for c in match if c.isdigit()])
    return [match,seq,operator,base_length,read_seq,insert_bases,items[0],read_length_overlapped_with_target_region,flag]

def extract_feature_matrix(reads_in_sam):
    '''extract match_maatrix, base_matrix and insert_matrix from sam file'''
    global base_matrix, match_matrix,insert_bases_matrix,reads_ID,reads_flag
    match_matrix,base_matrix,insert_bases_matrix,reads_ID,reads_flag  = ([] for i in range(5))
    for line in reads_in_sam:
        parsed_read = parse_sam(line)
        if parsed_read[7] > args.lengthFilter * len_ref:
            match_matrix.append(parsed_read[0])
            base_matrix.append(parsed_read[4])
            insert_bases_matrix.append(parsed_read[5])
            reads_ID.append(parsed_read[6])
            reads_flag.append(parsed_read[8])
    match_matrix = map(list,zip(*match_matrix))
    base_matrix = map(list,zip(*base_matrix))
    insert_bases_matrix = map(list,zip(*insert_bases_matrix))

def get_allele_freq_lst():
    global ct_base_lst, consensus_allele, depth, allele_freq_lst, ct_base_lst
    ct_base_lst,ct_base,allele_freq_lst,depth,consensus_allele,allele_freq = ([] for i in range(6))
    base_lst =  ["A","T","C","G","D"]
    for i in range(len(base_matrix)):
        ct_miscall=0
        ct_base = [0 for i in range(5)] 
        allele_freq = [0 for i in range(5)]
        for base in range(len(base_lst)):
            ct_base[base] = base_matrix[i].count(base_lst[base])
        total_ct = sum(ct_base)
        depth.append(total_ct)
        for base in range(len(base_lst)):
            allele_freq[base] = ct_base[base]/total_ct
        for freq in range(len(allele_freq)):
            if allele_freq[freq] < args.mismatch:
                ct_miscall += ct_base[freq]
        ct_base_lst.append(ct_base)
        allele_freq_lst.append(allele_freq)
        max_af = max(allele_freq)
        idx = allele_freq.index(max_af)
        consensus_allele.append('ATCGD'[idx])

def get_indel_rate():
    global del_rate, ins_rate,ct_del_lst, ct_ins_lst, ct_dep_lst
    ct_del_lst, ct_ins_lst,  ct_dep_lst,del_rate, ins_rate = ([] for i in range(5))
    for i in range(len(match_matrix)):
        match_count = match_matrix[i].count('1')
        del_count = match_matrix[i].count('0')
        ins_count = match_matrix[i].count('2')
        base_depth = match_count + del_count + ins_count
        del_rate.append(del_count/float(base_depth))
        ins_rate.append(ins_count/float(base_depth))
        ct_del_lst.append(del_count)
        ct_ins_lst.append(ins_count)
        ct_dep_lst.append(base_depth)

def get_consensus_insert_bases():
    global consensus_insert,snv_ins_freq
    consensus_insert,sninsert_freq,A_insert_freq,T_insert_freq,C_insert_freq,G_insert_freq =( [] for i in range(6))
    sninsert_freq = []
    for i in range(len(insert_bases_matrix)):
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
    snv_ins_freq  = [A_insert_freq,T_insert_freq,C_insert_freq,G_insert_freq]    
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

def large_deletion_finder():
    global true_large_del
    connected_deletion_matrix = [['0' for i in range(len(match_matrix[0]))] for j in range(len(match_matrix)-1)] 
    connected_deletion_freq =[0.0 for i in range(len(match_matrix)-1)]
    large_del_pos = []
    large_del_freq = []
    large_del = []
    true_large_del = []
    for i in range(len(match_matrix)-1):
        for j in range(len(match_matrix[0])):
            if match_matrix[i][j] == '0' and match_matrix[i+1][j] == '0':
                connected_deletion_matrix[i][j] = '1'
            if match_matrix[i][j] == '.' or match_matrix[i+1][j] == '.':
                connected_deletion_matrix[i][j] = '.'
    for i in range(len(match_matrix)-1):
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

def replace_confusing_deletion(i,j):
    '''correct deletion misplaced by aliner'''
    del_len = 1
    while match_matrix[i+del_len][j] == '0':
        del_len+=1
    if (i+2*del_len<len(base_matrix)
        and i+del_len not in true_large_del):
        for k in range(del_len):
            if base_matrix[i+k+del_len][j]!=consensus_allele[i+k]:
                break
            else:
                base_matrix[i+k][j] = consensus_allele[i+k]
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
        if consensus_allele[i] == 'D' and args.deletion + args.del_homo_p*(ref_context_info[i][2]-1) >=1:
            consensus_allele[i] = ref_seq[i]
        for j in range(len(base_matrix[i])):
            if match_matrix[i][j] == '0':        
                adjust_deletion_pos(i,j)
                if match_matrix[i][j] == '0':
                    replace_confusing_deletion(i,j)
                    if match_matrix[i][j] == '0':
                        base_matrix[i][j] = consensus_allele[i]
                        match_matrix[i][j] = '1'
    elif (i>0 
          and del_rate[i] >=  args.deletion + args.del_homo_p*(ref_context_info[i][2]-1)
          and del_rate[i] <= args.deletion + args.del_homo_p*ref_context_info[i][2]):
        for j in range(len(base_matrix[i])):
            if match_matrix[i][j] == '0':
                correct_del_in_new_2bp_homopolymer(i,j)
    elif del_rate[i] >= 1 - args.insert:
        for j in range(len(base_matrix[i])):
            base_matrix[i][j] = 'D'    
            match_matrix[i][j] = '0'

def correct_single_insert(i):
    for ins in ['A','T','C','G']:
        if snv_ins_freq['ATCG'.index(ins)][i]<= args.insert +args.ins_homo_p*(ref_context_info[i][2]-1):
            for j in range(len(insert_bases_matrix[i])):
                if insert_bases_matrix[i][j] == ins:
                    if (i>0 
                        and base_matrix[i-1][j] 
                        and allele_freq_lst[i-1]['ATCGD'.index(base_matrix[i-1][j])] <= args.mismatch 
                        and  insert_bases_matrix[i][j] == consensus_allele[i-1]):
                        base_matrix[i-1][j]=consensus_allele[i-1]
                        match_matrix[i-1][j] = '1'
                    elif (i<len(base_matrix)-1 
                          and base_matrix[i+1][j] 
                          and allele_freq_lst[i+1]['ATCGD'.index(base_matrix[i+1][j])] <= args.mismatch 
                          and  insert_bases_matrix[i][j] == consensus_allele[i+1]):
                        base_matrix[i+1][j]=consensus_allele[i+1]
                        match_matrix[i+1][j] = '1'
                    insert_bases_matrix[i][j] = 'O'
                    match_matrix[i][j] = '1'
   
def correct_large_insert(i):
    if ins_rate[i] <= args.insert +args.ins_homo_p*(ref_context_info[i][2]-1):
        for j in range(len(insert_bases_matrix[i])):
            if insert_bases_matrix[i][j] not in ['','O']:
                if (i>0 
                   and base_matrix[i-1][j] not in ['.','R'] 
                   and allele_freq_lst[i-1]['ATCGD'.index(base_matrix[i-1][j])] <= args.mismatch 
                   and  insert_bases_matrix[i][j] == consensus_allele[i-1]):
                    base_matrix[i-1][j]=consensus_allele[i-1]
                    match_matrix[i-1][j] = '1'
                elif (i<len(base_matrix)-1 
                      and base_matrix[i+1][j] not in ['.','R'] 
                      and allele_freq_lst[i+1]['ATCGD'.index(base_matrix[i+1][j])] <= args.mismatch 
                      and  insert_bases_matrix[i][j] == consensus_allele[i+1]):
                    base_matrix[i+1][j]=consensus_allele[i+1]
                    match_matrix[i+1][j] = '1'
                elif len(insert_bases_matrix[i][j]) == 1:
                    if (base_matrix[i][j] not in ['.','R'] 
                       and allele_freq_lst[i]['ATCGD'.index(insert_bases_matrix[i][j])] >= args.mismatch):
                        base_matrix[i][j]= insert_bases_matrix[i][j]
                insert_bases_matrix[i][j] = 'O'
                match_matrix[i][j] = '1'
    elif ins_rate[i] >= 1 - args.deletion:
        for j in range(len(insert_bases_matrix[i])):
            if insert_bases_matrix[i][j] == 'O':
                insert_bases_matrix[i][j] =consensus_insert[i]
                match_matrix[i][j] = '2'
  
def correct_miscall():
    for i in range(len(allele_freq_lst)):
        for j in range(len(base_matrix[i])):
            if  (base_matrix[i][j] not in ['.','R','D'] 
                 and (allele_freq_lst[i]['ATCGD'.index(base_matrix[i][j])] <= args.mismatch)):
                if (i>0 and base_matrix[i-1][j] not in ['.','R','D']
                    and allele_freq_lst[i-1]['ATCG'.index(base_matrix[i][j])] >= args.mismatch):
                    base_matrix[i-1][j]=base_matrix[i][j]
                elif (i<len(base_matrix)-1
                    and base_matrix[i+1][j] not in ['.','R','D']
                    and allele_freq_lst[i+1]['ATCG'.index(base_matrix[i][j])] >= args.mismatch):
                    base_matrix[i+1][j]=base_matrix[i][j]
                base_matrix[i][j] = consensus_allele[i]
                match_matrix[i][j] = '1'

def get_reads_from_input(input_file):
    reads_in_sam=[]
    with open (input_file,"r") as sam_file:
        for line in sam_file.readlines():
           if line[0] != '@' and random.randrange(1,101)/100.0 <= args.downsample_freq:
               line1 = line.split('\t')
               if line1[4] != '0' and str(line1[2]) == region_chr and line1[1] in ['0','16','2048','2064']: 
                  reads_in_sam.append(line)
    return reads_in_sam

def correct_reads(reads_in_sam):
    global mis_rate, allele_freq_lst, snv_ins_freq,base_matrix, match_matrix,insert_bases_matrix,reads_ID,reads_flag,mis_rate0,ins_rate0,del_rate0,ct_base_lst0,ct_ins_lst0,ct_dep_lst0 ,ct_del_lst0
    correction,allele_freq_lst,base_matrix, match_matrix,insert_bases_matrix,reads_ID,reads_flag,mis_rate  = ([] for i in range(8))
    extract_feature_matrix(reads_in_sam)
    get_allele_freq_lst()
    get_consensus_insert_bases()
    for i in range(len(allele_freq_lst)):
        mis_rate.append(1-allele_freq_lst[i]["ATCG".index(ref_seq[i])]-allele_freq_lst[i][4])
    get_indel_rate() 
    large_deletion_finder()
    mis_rate0 = mis_rate
    del_rate0 = del_rate
    ins_rate0 = ins_rate
    ct_base_lst0 = ct_base_lst
    ct_ins_lst0 = ct_ins_lst
    ct_del_lst0 = ct_del_lst
    ct_dep_lst0 = ct_dep_lst
    for i in range(len(allele_freq_lst)):
        correct_large_insert(i)
        correct_deletion(i)
    get_indel_rate()
    for i in range(len(allele_freq_lst)):
        correct_single_insert(i)
    get_indel_rate()
    get_allele_freq_lst()
    correct_miscall()
    match_matrix = map(list,zip(*match_matrix)) 
    base_matrix = map(list,zip(*base_matrix))
    insert_bases_matrix = map(list,zip(*insert_bases_matrix))
    for i in range(len(match_matrix)):
        corrected_seq = ''
        for j in range(len(match_matrix[i])):
            if match_matrix[i][j] == '1':
                corrected_seq += base_matrix[i][j]
            elif match_matrix[i][j] == '2':
                corrected_seq += (insert_bases_matrix[i][j]+base_matrix[i][j])
        corrected_seq = corrected_seq.replace('D','')
        corrected_seq = corrected_seq.replace('R','')
        if reads_flag[i] == '0' or reads_flag[i] == '2048':
            correction.append(corrected_seq)
        else:
            correction.append(RevC(corrected_seq))
    return [reads_ID,correction]

def output_corrected_reads(cor_reads, input_file):
    with open(input_file[:-3]+"corrected.fasta","w+") as output_file:
        for i in range(len(cor_reads[1])):
            output_file.write('>' + cor_reads[0][i] + '\n' + cor_reads[1][i] + '\n')

def output_error_rate(input_file):
    with open(input_file[:-3]+"error_rate_file_before","w+") as error_rate_file:
        for i in range(len(mis_rate)):
            error_rate_file.write('\t'.join(map(str,[region_chr, i + region_start, round(mis_rate0[i],5),round(del_rate0[i],5),round(ins_rate0[i],5)]+map(str,ct_base_lst0[i])+map(str,[ct_del_lst0[i],ct_ins_lst0[i],ct_dep_lst0[i]])+ ref_context_info[i]))+'\n')

    with open(input_file[:-3]+"error_rate_file_after","w+") as error_rate_file:
        for i in range(len(mis_rate)):
            error_rate_file.write('\t'.join(map(str,[region_chr, i + region_start, round(mis_rate[i],5),round(del_rate[i],5),round(ins_rate[i],5)]+map(str,ct_base_lst[i])+map(str,[ct_del_lst[i],ct_ins_lst[i],ct_dep_lst[i]])+ ref_context_info[i]))+'\n')
            
def main():
    global region_chr, region_start, region_end,args 
    parser = get_argument()
    args = parser.parse_args()
    samfile = pysam.AlignmentFile(args.input,"r")
    region_chr = args.targetRegion[:args.targetRegion.find(":")]
    region_start = int(args.targetRegion[args.targetRegion.find(":")+1:args.targetRegion.find("-")])
    region_end = int(args.targetRegion[args.targetRegion.find("-")+1:])
    print "#################################################################################"
    print "Correcting reads in region" + args.targetRegion
    parse_fasta(args.reference)
    time_b = time.time()
    print ("ALEC took %s seconds to get reference info" % ( round((time_b - start_time),0)))
    reads = get_reads_from_input(args.input)
    time_c = time.time()
    print ("ALEC took %s seconds to get reads" % ( round((time_c - time_b),0)))
    cor_reads = correct_reads(reads)
    output_error_rate(args.input)
    output_corrected_reads(cor_reads, args.input)
    print ("ALEC took %s seconds to correct" % (round((time.time() - time_c),0)))
    print ("ALEC took %s seconds to correct %s reads in region %s" % ( round((time.time() - start_time),2), len(cor_reads[1]), args.targetRegion))
    print "#################################################################################"

if __name__ == "__main__":
    main()


