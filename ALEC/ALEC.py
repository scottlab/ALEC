#!/usr/bin/env python
#title           :pyscript.py
#description     :correct long sequencing(pacbio) reads after aglinment
#author          :Yao Yang
#date            :20150530
#version         :0.1
#usage           :python alec.py
#notes           :
#python_version  :2.7.2
#==============================================================================

'''Import the modules needed to run the script'''
import sys
import re
import time
import argparse
'''Define global variables'''
start_time = time.time()
num_del = 0
num_ins = 0
num_mis = 0
total_num = 0
ref_seq = ''
len_ref = 0

'''Open input and output files'''
try:
    input_file = open(sys.argv[1],"r")
    reference_file = open(sys.argv[2],"r")
    output_file = open(sys.argv[1][:-3]+'corrected.fasta',"w+")
except:
    print("""
ALEC takes a fasta file as reference and a SAM file as the raw data alignment information, and automatically generates a corrected fasta file as output in the same directory as raw data.

The usage of the script is as below:

    python ALEC.py input.sam reference.fasta

Note: the script only takes one single sequence as reference each time.
    """)
    sys.exit(0)

'''define function for getting DNA sequence from reference fasta file'''
def parse_fasta(reference_file):
    global ref_seq
    counter = 1
    ref_context_info = []
    for line in reference_file:
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
    return [ref_seq,ref_context_info]

'''define function for getting DNA reverse complement'''
def RevC(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N','D':'D'}
    return "".join([seq_dict[base] for base in reversed(seq)])

'''define funcion for extracting single read information from each line in sam file'''
def parse_sam(read):
    global total_num, ref_seq
    len_ref = len(ref_seq)
    read = read.rstrip()
    items = read.split('\t')
    [start,cigar,seq] = [int(items[3]),items[5],items[9]]
    match = '.'*(int(items[3]) - 1)
    insert_bases = ['']*(int(items[3]) - 1)
    num,align = ([] for i in range(2))
    insert = False
    read_seq = '.'*(int(items[3]) - 1)
    while len(cigar)>0:    
        for i in range(len(cigar)):
            if cigar[i].isalpha():
                pos = i
                break
        num.append(int(cigar[:pos]))
        align.append(cigar[pos])
        cigar = cigar[pos+1:]
    if align[0] in ['S','H']:
        seq = seq[num[0]:]
        align.pop(0)
        num.pop(0)
    if align[-1] in ['S','H']:
        seq = seq[:len(seq) - num[-1]]
        align.pop()
        num.pop()
    seq_cy = seq
    for i in range(len(num)):
        if align[i] == 'M':
            read_seq += seq_cy[:num[i]]
            seq_cy = seq_cy[num[i]:]
            if not insert:
                match += '1'*(num[i])
                insert_bases += ['O']*(num[i]) 
            else:
                match += '1'*(num[i]-1)
                insert_bases += ['O']*(num[i] - 1)
                insert = False
        elif align[i] == 'D':
            read_seq += 'D'*num[i]
            match += '0'*num[i]
            insert_bases += ['']*(num[i])
        elif align[i] == 'I':
            match += '2'
            insert_bases.append(seq_cy[:num[i]])
            seq_cy = seq_cy[num[i]:]
            insert = True 
    total_num += (len(match) - (int(items[3]) - 1))
    match += '.'*(len_ref - len(match))
    read_seq += (len_ref - len(read_seq))*'.'
    insert_bases += [''] * (len_ref - len(insert_bases))
    return [match,seq,align,num,read_seq,insert_bases,items[0]]

'''define function for extract match_maatrix, base_matrix and insert_matrix from sam file'''
def extract_feature_matrix(reads_in_sam):
    match_matrix,base_matrix,insert_bases_matrix,reads_ID  = ([] for i in range(4))
    for line in reads_in_sam:
        parse_sam_output = parse_sam(line)
        match_matrix.append(parse_sam_output[0])
        base_matrix.append(parse_sam_output[4])
        insert_bases_matrix.append(parse_sam_output[5])
        reads_ID.append(parse_sam_output[6])
    match_matrix = map(list,zip(*match_matrix))
    base_matrix = map(list,zip(*base_matrix))
    insert_bases_matrix = map(list,zip(*insert_bases_matrix))
    return [match_matrix, base_matrix, insert_bases_matrix,reads_ID]

def get_allele_freq_lst(base_matrix):
    allele_freq_lst = []
    most_possible_allele = ''
    for i in range(len(base_matrix)):
         A_ct = base_matrix[i].count('A')
         T_ct = base_matrix[i].count('T')
         C_ct = base_matrix[i].count('C')
         G_ct = base_matrix[i].count('G')
         D_ct = base_matrix[i].count('D')
         total_ct = float(A_ct + T_ct + C_ct + G_ct + D_ct)
         allele_freq = [A_ct/total_ct,T_ct/total_ct,C_ct/total_ct,G_ct/total_ct,D_ct/total_ct]
         allele_freq_lst.append(allele_freq)
         max_af = max(allele_freq)
         idx = allele_freq.index(max_af)
         most_possible_allele += 'ATCGD'[idx]
    return [allele_freq_lst,most_possible_allele]

def get_indel_rate(match_matrix):
    del_rate = []
    ins_rate = []
    for i in range(len(match_matrix)):
        match_count = match_matrix[i].count('1')
        del_count = match_matrix[i].count('0')
        ins_count = match_matrix[i].count('2')
        total_reads = match_count + del_count + ins_count
        del_rate.append(del_count/float(total_reads))
        ins_rate.append(ins_count/float(total_reads))
    return [del_rate,ins_rate]

def get_most_possbile_insert_bases(insert_bases_matrix):
    most_possible_insert = []
    for i in insert_bases_matrix:
        unique_insert =list(set(i))
        if '' in unique_insert:
            unique_insert.remove('')
        ct_insert_bases = []
        for element in unique_insert:
            ct_insert_bases.append(i.count(element))
        if ct_insert_bases == []:
            most_possible_insert.append('')
        elif max(ct_insert_bases)/float(sum(ct_insert_bases)) >= 0.8:
            if unique_insert[ct_insert_bases.index(max(ct_insert_bases))] != 'O':
                most_possible_insert.append(unique_insert[ct_insert_bases.index(max(ct_insert_bases))])
            else:
                most_possible_insert.append('')
        else:
            most_possible_insert.append('')
    return most_possible_insert

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

def correct_deletion(base_matrix, del_rate, most_possible_allele,match_matrix):
    global num_del
    for i in range(len(del_rate)):
        if del_rate[i] == 0:
            pass
        elif del_rate[i] <= 0.15 or del_rate[i] <= 0.15+0.05*(ref_context_info[i][2]-2):
            for j in range(len(base_matrix[i])):
                if base_matrix[i][j] == 'D':
                    base_matrix[i][j] = most_possible_allele[i]
                    match_matrix[i][j] = '1'
                    num_del += 1
    return [base_matrix,match_matrix]

def correct_insert(base_matrix, ins_rate, insert_bases_matrix,most_possible_allele,match_matrix,allele_freq_lst):
    global num_ins
    polybase = []
    for i in range(len(allele_freq_lst)):
        if sorted(allele_freq_lst[i])[-2]>=0.15:
            polybase.append(i)
    for i in range(len(ins_rate)):
        if ins_rate[i] == 0:
            pass
        elif i in polybase:
            for j in range(len(insert_bases_matrix[i])):
                if insert_bases_matrix[i][j] == 'ATCGD'[allele_freq_lst[i].index(sorted(allele_freq_lst[i])[-2])] or  insert_bases_matrix[i][j] == most_possible_allele[i]:
                    match_matrix[i][j] = '1'
                    base_matrix[i][j] =  insert_bases_matrix[i][j]
                    insert_bases_matrix[i][j] = 'O'
                    num_ins += 1
                elif insert_bases_matrix[i-1][j] == base_matrix[i-1][j] == 'ATCGD'[allele_freq_lst[i].index(sorted(allele_freq_lst[i])[-2])] or  insert_bases_matrix[i-1][j] == base_matrix[i-1][j] == most_possible_allele[i]:
                    match_matrix[i-1][j] = '1'
                    base_matrix[i][j] =  insert_bases_matrix[i-1][j]
                    insert_bases_matrix[i-1][j] = 'O'
                    num_ins += 1
                elif insert_bases_matrix[i+1][j] == 'ATCGD'[allele_freq_lst[i].index(sorted(allele_freq_lst[i])[-2])] or  insert_bases_matrix[i+1][j] == most_possible_allele[i]:
                    match_matrix[i+1][j] = '1'
                    base_matrix[i][j] =  insert_bases_matrix[i+1][j]
                    insert_bases_matrix[i+1][j] = 'O'
                    num_ins += 1
        elif ins_rate[i] <= 0.15 or ins_rate[i] <= 0.15+0.05*(ref_context_info[i][2]-2):
            for j in range(len(insert_bases_matrix[i])):
                if insert_bases_matrix[i][j] not in ['','O']:
                    insert_bases_matrix[i][j] = 'O'
                    match_matrix[i][j] = '1'
                    num_ins += 1
    return [base_matrix,insert_bases_matrix,match_matrix]

def correct_base(match_matrix,most_possible_allele,base_matrix,allele_freq_lst,ins_rate,most_possible_insert,insert_bases_matrix):
    global num_mis, num_ins
    for i in range(len(allele_freq_lst)):
        if allele_freq_lst[i][4] >= 0.75:
            for j in range(len(base_matrix[i])):
                if base_matrix[i][j] != 'D':
                    base_matrix[i][j] = 'D'
                    match_matrix[i][j] = '0'
                    num_ins += 1
        else:
            for j in range(len(base_matrix[i])):
                if  (base_matrix[i][j] not in ['.','R']) and (allele_freq_lst[i]['ATCGD'.index(base_matrix[i][j])] <= 0.15):
                    base_matrix[i][j] = most_possible_allele[i]
                    num_mis += 1
        if ins_rate[i] >= 0.75:
            for j in range(len(insert_bases_matrix[i])):
                if insert_bases_matrix[i][j] == 'O':
                    insert_bases_matrix[i][j] =most_possible_insert[i]
                    match_matrix[i][j] = '2'
                    num_del += 1
    return [base_matrix,match_matrix,insert_bases_matrix]  

def correct_sam(sam_file):
    global len_ref, ref_context_info
    mis_rate = []
    reads_in_sam,seq,strand,cigar,start_pos,operation,length,correction,allele_freq_lst  = ([] for i in range(9))
    for line in sam_file:
        if line[0] != '@':
            line1 = line.split('\t')
            if line1[4] != '0' and line1[1] in ['0','16']:
                reads_in_sam.append(line)
    for read in reads_in_sam:
        read1 = read.split('\t')
        parse_sam_output = parse_sam(read)
        seq.append(parse_sam_output[1])
        operation.append(parse_sam_output[2])
        length.append(parse_sam_output[3])
        strand.append(read1[1])
        cigar.append(read1[5])
        start_pos.append(int(read1[3]))
    [reference,ref_context_info] = parse_fasta(reference_file)
    len_ref = len(reference)
    match_matrix, base_matrix, insert_bases_matrix, reads_ID = extract_feature_matrix(reads_in_sam)
    most_possible_insert = get_most_possbile_insert_bases(insert_bases_matrix)
    allele_freq_lst,most_possible_allele = get_allele_freq_lst(base_matrix)
    for i in range(len(allele_freq_lst)):
        mis_rate.append(1-allele_freq_lst[i]["ATCG".index(reference[i])]-allele_freq_lst[i][4])
    del_rate,ins_rate = get_indel_rate(match_matrix) 
    [base_matrix,match_atrix] = correct_deletion(base_matrix, del_rate, most_possible_allele, match_matrix)
    [base_matrix,insert_bases_matrix,match_matrix] = correct_insert(base_matrix, ins_rate, insert_bases_matrix,most_possible_allele,match_matrix,allele_freq_lst)
    allele_freq_lst,most_possible_allele = get_allele_freq_lst(base_matrix)
    [base_matrix,match_matrix,insert_bases_matrix] = correct_base(match_matrix,most_possible_allele,base_matrix,allele_freq_lst,ins_rate,most_possible_insert,insert_bases_matrix)
    
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
        if strand[i] == '0':
            correction.append(corrected_seq)
        else:
            correction.append(RevC(corrected_seq))
    return [reads_ID,correction]

cor_reads = correct_sam(input_file)

print (">>>>>>>> %s seconds <<<<<<<<" % (time.time() - start_time))
print ("deletion rate: %s" %  (num_del/float(total_num)))
print ("insert rate: %s" % (num_ins/float(total_num)))
print ("mismatch rate: %s" % (num_mis/float(total_num)))
for i in range(len(cor_reads[1])):
    output_file.write('>' + cor_reads[0][i] + '\n' + cor_reads[1][i] + '\n')
    
