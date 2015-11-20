

# Import                                                                                    
import os
from NeedlemanWunsch import *
import random
import operator

def initial_pairwise_alignments(seqs, num_seq, subst, gap):
    # Create initial population from pairwise sequence alignments                               
    Alignments = {}
    n = 0
    for i in range(num_seq-1):
        for j in range(i+1, num_seq):
            seq1 = seqs[i]
            seq2 = seqs[j]
            (aligned1, aligned2) = nw(seq1, seq2, subst, gap)
            Alignments[(i, n)] = aligned1
            Alignments[(j, n)] = aligned2
            n += 1
    return Alignments


# Calculate fitness based on sum of pairwise alignment scores                               
def calc_fitness(organism, subst, gap, max_len):
    total_score = 0
    for i in range(len(organism)):
        for j in range(i+1, len(organism)):
            score = 0
            for k in range(max_len):
                if operator.xor(bool(organism[i][k] == '-'), bool(organism[j][k] == '-')):
                    if organism[i][k] != organism[i][j]:
                        score += gap
                elif organism[i][k] != '-' and organism[j][k] != '-':
                    score += subst[organism[i][k], organism[j][k]]
            total_score += score
    return total_score

def initial_population(Alignments, pop_size, num_seq, max_offset, max_len):
    Population = {}
    # People a population & score                                                               
    for i in range(pop_size):
        for j in range(num_seq):
            rand_index = random.randint(0, num_seq-1)
            while not (j, rand_index) in Alignments.keys():
                rand_index = random.randint(0, num_seq-1)
            chosen_seq = Alignments[(j, rand_index)]
            chosen_seq_l = list(chosen_seq)
            offset = random.randint(0, max_offset)
            for k in range(offset):
                chosen_seq_l.insert(0, '-')
            while len(chosen_seq_l) < max_len:
                chosen_seq_l.append('-')
            if i in Population.keys():
                Population[i].append(''.join(chosen_seq_l))
            else:
                Population[i] = [''.join(chosen_seq_l)]
    return Population


def initial_fitness(Population, pop_size, subst, gap, max_len):
    Fitness = {}
    for i in range(pop_size):
        Fitness[i] = calc_fitness(Population[i], subst, gap, max_len)
    return Fitness

def horizontal_recombination(D, mom, dad, num_seq):
    offspring = []
    for i in range(num_seq):
        P_prob = random.randint(0, 1)
        if P_prob == 0:
            offspring.append(D[mom][i])
        else:
            offspring.append(D[dad][i])
    return offspring

def vertical_recombination(D, mom, dad, l, num_seq, max_len):
    offspring = []
    for i in range(num_seq):
        m_offset = 0
        d_offset = 0
        while D[mom][i][m_offset] == '-':
            m_offset += 1
        d_offset = 0
        while D[dad][i][d_offset] == '-':
            d_offset += 1
        mom_string = ''.join(D[mom][i])
        dad_string = ''.join(D[dad][i])
        for j in range(m_offset, m_offset+ l):
            if D[mom][i][j] == '-':
                m_offset += 1
        for k in range(d_offset, d_offset+l):
            if D[dad][i][k] == '-':
                d_offset += 1
        seq = mom_string[:l+m_offset] + dad_string[l+d_offset:]
        seq_list = list(seq)
        while len(seq_list) < max_len:
            seq_list.append('-')
        offspring.append(''.join(seq_list))
    return offspring

def gap_extension(D, org, num_seq, max_len):
    gap_blocks = find_gap_blocks(D, org, num_seq, max_len)
    new_org = []
    block = random.randint(0, len(gap_blocks)-1)
    pos = random.randint(gap_blocks[block][0],gap_blocks[block][0]+ gap_blocks[block][1])
    for i in range(num_seq):
        M_seq_l = list(D[org][i])
        M_seq_l.insert(pos, '-')
        new_org.append(''.join(M_seq_l))
    return new_org

def gap_addition(D, org, num_seq, max_len):
    new_org = []
    pos = random.randint(0, max_len-1)
    for i in range(num_seq):
        seq_l = list(D[org][i])
        seq_l.insert(pos, '-')
        new_org.append(''.join(seq_l))
    return new_org

def gap_reduction(D, org, num_seq, max_len):
    probs = []
    new_org = []
    blocks = find_gap_blocks(D, org, num_seq, max_len)
    if len(blocks) == 0:
        return D[org]
    for i in range(len(blocks)):
        for j in range(int(max_len / blocks[i][1])):
            probs.append(blocks[i][0])
    pos = random.choice(probs)
    for j in range(num_seq):
        seq_l = list(D[org][j])
        del seq_l[pos]
        seq_l.append('-')
        new_org.append(''.join(seq_l))
    return new_org

def find_gap_blocks(D, org, num_seq, max_len):
    gaps = []
    blocks = []
    for i in range(max_len):
        n = []
        for j in range(num_seq):
            n.append(D[org][j][i])
        if n[1:] == n[:-1] and n[0] == '-':
            gaps.append(i)
    block_len = 1
    block_start = gaps[0]
    for k in range(len(gaps)):
        if k+1 < len(gaps) and gaps[k]+1 == gaps[k+1]:
            block_len += 1
        else:
            blocks.append((block_start, block_len))
            if k+1 < len(gaps):
                block_start = gaps[k+1]
            block_len = 1
    return blocks

def trim_gaps(D, org, num_seq, max_len):
    blocks = find_gap_blocks(D, org, num_seq, max_len)
    new_org = []
    for i in range(num_seq):                                        
        seq = []
        # Blocks at front and back
        if blocks[0][0] == 0:
            no_gap_start = blocks[0][0]+blocks[0][1]
            no_gap_end = blocks[-1][0]
        # Blocks only at back
        else:
            no_gap_start = 0
            no_gap_end = blocks[-1][0]
        seq_s = D[org][i][no_gap_start:no_gap_end]
        seq.append(seq_s)
        new_org.append(''.join(seq[0]))
    return new_org

