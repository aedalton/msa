# Implement genetic algorithm for multple sequence alignment                                
# Inspired by paper "A simple genetic algorithm for multiple sequence alignment"            
# (C. Gondro & B.P. Kinghorn)                                                               

# Import                                                                                    
import os
from NeedlemanWunsch import *
from HelperFunctions import *
import random
import operator

# For now, just take first three sequences from globin_fragments.fasta, use defined         
# substitution table, gap penalty, etc.                                                     
# Later, get input from user                                                                
# seqs: the path to a FASTA file that contains the sequences to be aligned
seqs = FASTAFile('globin_fragments.fasta')
# subst: The substitution matrix to be used to align
# For amino acid sequences, can use blosum62, blosum45, or an exact 20x20 matrix
# For nucleotide sequences, can use an exact 4x4 matrix
subst = blosum62
# gap: Gap penalty for inserting any gap in the sequence (not a linear model)
gap = -4                                                  
# num_seq: Number of sequences to be aligned (first n in the file provided)
num_seq = 3
# pop_size: Size of the population
pop_size = 25
# generations: Number of generations
generations = 500
# culling_percentage: Bottom percentage to be removed from the population without 
# reproducing for every generation
culling_percentage = 0.4
# num_least_fit: The number of organisms in the population that are removed
# without reproducing
num_least_fit = int(culling_percentage*pop_size)
# HR_prob: The probability of a horizontal recombination event occurring
HR_prob = 0.5
# VR_prob: The probability of a vertical recombination event occurring
VR_prob = 0.3
# GE_prob: The probability of a gap extension mutation occurring
GE_prob = 0.05
# GA_prob: The probability of a gap addition mutation occurring
GA_prob = 0.1
# GR_prob: The probability of a gap reduction mutation occurring
GR_prob = 0.05


# Dictionary with tuple as key to hold all possible alignments                              
Alignments = initial_pairwise_alignments(seqs, num_seq, subst, gap)

# Calculate properties of sequences in Alignments
# max_seq_len: The length of the longest of the sequences in the initial
# population
max_seq_len = 0
for i in range(num_seq):
    if len(seqs[i]) > max_seq_len:
        max_seq_len = len(seqs[i])
# max_offset: The maximum number of gaps to be inserted at the beginning of the
# sequence, defined to be 20% of the length of the longest ungapped sequence
max_offset = int(0.2*max_seq_len)
# max_len: The maximum length of the sequences, defined to be 135% of the length
# of the longest ungapped sequence, plus the calculated offset
max_len = int(max_seq_len*1.35) + max_offset
# VR_index: To be used in vertical recombination events, the point at which the
# sequences are broken and conjoined 
VR_index = random.randint(0, max_seq_len)

# Dictionary with a tuple as a key to hold all organisms                                    
Population = initial_population(Alignments, pop_size, num_seq, max_offset, max_len)

# Dictionary with an integer as a key to the fitnes score of the coressponding              
# organism in the population                                                                
Fitness = initial_fitness(Population, pop_size, subst, gap, max_len)

# Generations
for gen in range(generations):
    # Survival of the fitest
    # Cull the bottom percentage of the population before reproduction events
    for i in range(num_least_fit):
        del Population[sorted(Fitness.items(), key=operator.itemgetter(1), reverse = True)[-1][0]]
        del Fitness[sorted(Fitness.items(), key=operator.itemgetter(1), reverse = True)[-1][0]]
    # Horizontal & Vertical Recombination
    # Create the same number of offspring that was just removed from the population
    # to maintain population size
    for i in range(num_least_fit):
    # Select 2 unique parents at random to create an offspring                                     
        mom = random.choice(Population.keys())
        dad = random.choice(Population.keys())
        while mom == dad:
            dad = random.choice(Population.keys())
        # Choose a random floating point number between 0 and 1 to determine what type
        # of reproduction (horizontal recombimation, vertical recombination, or copy)
        # will occur
        R_prob = random.random()
        # Horizontal recombination
        if R_prob < HR_prob:
            offspring = horizontal_recombination(Population, mom, dad, num_seq)
        # Vertical recombination
        elif R_prob > HR_prob and R_prob < HR_prob+VR_prob:
            offspring = vertical_recombination(Population, mom, dad, VR_index, num_seq, max_len)
        # Copy from mom or dad
        else:
            P_prob = random.randint(0, 1)
            if P_prob == 0:
                offspring = Population[mom]
            else:
                offspring = Population[dad]
        # Insert with a unique organism index into the population and calculate fitness
        key = sorted(Population)[-1] + 1
        Population[key] = offspring
        Fitness[key] = calc_fitness(Population[key], subst, gap, max_len)

    # Mutation
    # Any organism has the same small probability of being mutated
    for i in range(pop_size):
        # Select organism to be mutated
        old_org = random.choice(Population.keys())
        # Choose a random floating point number between 0 and 1 to determine what type
        # of mutation may occur (gap extension, gap addition, or gap reduction)
        M_prob = random.random()
        # Gap Extension
        if M_prob < GE_prob:
            # Extend gap
            new_org = gap_extension(Population, old_org, num_seq, max_len)
            # Remove old organism from population
            del Population[old_org]
            del Fitness[old_org]
            # Insert mutated new organism into population and calculate its fitness
            key = sorted(Population)[-1] + 1
            Population[key] = new_org
            Fitness[key] = calc_fitness(Population[key], subst, gap, max_len)
            # Update maximum length and adjust length of other organisms in the population
            max_len += 1
            for k in Population:
                if k != key:
                    for j in range(num_seq):
                        l = list(Population[k][j])
                        l.append('-')
                        s = ''.join(l)
                        Population[k][j] = s

        # Gap Addition
        elif M_prob > GE_prob and M_prob < GE_prob + GA_prob:
            # Insert gap
            new_org = gap_addition(Population, old_org, num_seq, max_len)
            # Remove old organism from population
            del Population[old_org]
            del Fitness[old_org]
            # Insert mutated new organism into population and calculate its fitness
            key = sorted(Population)[-1] + 1
            Population[key] = new_org
            Fitness[key] = calc_fitness(Population[key], subst, gap, max_len)
            # Update maximum length and adjest length of other organisms in the population
            max_len += 1
            for k in Population:
                if k != key:
                    for j in range(num_seq):
                        l = list(Population[k][j])
                        l.append('-')
                        s = ''.join(l)
                        Population[k][j] = s

        # Gap Reduction
        elif M_prob > GE_prob + GA_prob and M_prob < GE_prob + GA_prob + GR_prob:
            # Reduce gap
            new_org = gap_reduction(Population, old_org, num_seq, max_len)
            # Remove old organism from population
            del Population[old_org]
            del Fitness[old_org]
            # Insert new mutated organism into population and calculate its fitness
            key = sorted(Population)[-1] + 1
            Population[key] = new_org
            Fitness[key] = calc_fitness(Population[key], subst, gap, max_len)
    # Report percent of evolution complete
    if gen % int(0.1 * generations) == 0:
        print str(100*gen/generations) + "%"
    
# Report most fit organism
most_fit = sorted(Fitness.items(), key=operator.itemgetter(1))[-1][0]
trimmed_most_fit = trim_gaps(Population, most_fit, num_seq, max_len)
final_fitness = Fitness[most_fit]
print 
print "MULTIPLE SEQUENCE ALIGNMENT"
for i in range(len(trimmed_most_fit)):
    print trimmed_most_fit[i]
print "SCORE"
print final_fitness

# Cleanup 
os.system("rm *c")
