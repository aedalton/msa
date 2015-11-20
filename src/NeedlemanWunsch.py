

# File contains functions to perform Needleman-Wunsch global sequence alignment

# Import                                                                                               
import numpy
from FASTAReader import *
from SubstitutionMatrix import *
from sys import argv, exit

# Define a class of objects that will hold score and backpointer for each entry in the matrix          
class Entry:
    def __init__(self, score, pointer):
        self.score = score
        self.pointer = pointer

# Function to create scoring matrix                                                                    
# Input: a and b, sequences to be aligned; subst, the substitution matrix to be used; gap, 
# the gap penalty to be used                           
# Output: returns scoring matrix                                                                       
def matrix(a, b, subst, gap):
    # Create an empty matrix of objects, m+1 by n+1                                                    
    M = numpy.empty([1+len(a), 1+len(b)], dtype='object')
    # Fill in matrix                                                                                   
    for i in range(len(a)+1):
        for j in range(len(b)+1):
            # Fill in first row and column                                                             
            if i==0 and j==0:
                # Upperleftmost entry does not have backpointer                                        
                M[i][j] = Entry(0, 'n')
            # First row                                                                                
            elif i==0:
                # First row is created by moving horizontally, levying a gap penalty for every move    
                s = M[i][j-1].score + gap
                M[i][j] = Entry(s, 'h')
            # First column                                                                             
            elif j==0:
                # First column is created by moving vertically, levying a gap penalty for every move   
                s = M[i-1][j].score + gap
                M[i][j] = Entry(s, 'v')
            # Fill in the remainder of the matrix                                                      
            else:
                # Calculate three possible scores (diagonal, vertical, horizontal)                     
                diagonal_score = M[i-1][j-1].score + subst[a[i-1], b[j-1]]
                vertical_score = M[i-1][j].score + gap
                horizontal_score = M[i][j-1].score + gap
                # Store highest score and pointer, prefering diagonal, then vertical, then horizontal  
                # Record maximum score with appropriate backpointer                              
                if max(diagonal_score, vertical_score, horizontal_score)==diagonal_score:
                    M[i][j] = Entry(diagonal_score, 'd')
                elif max(diagonal_score, vertical_score, horizontal_score)==vertical_score:
                    M[i][j] = Entry(vertical_score, 'v')
                elif max(diagonal_score, vertical_score, horizontal_score)==horizontal_score:
                    M[i][j] = Entry(horizontal_score, 'h')
    # Return matrix                                                                                    
    return M

# Function to backtrace through a scoring matrix to construct the best alignment                       
# Input: M, scoring matrix; a and b, sequences          
# Output: Best alignment string                                                                        
def backtrace (M, a, b):
    # Initialize lists for sequences                                                                   
    seq_1=[]
    seq_2=[]
    # Start at end of sequences                                                                        
    i = len(a)-1
    j = len(b)-1
    # While there is still a backpointer to follow (after all residues are alligned in case of global,
    # after score drops to zero for local alignments)                                                  
    while M[i][j].pointer != 'n':
        # If pointer is diagonal, "consume" one residue from each sequence                             
        if M[i][j].pointer == 'd':
            seq_1.append(a[i])
            seq_2.append(b[j])
            i -= 1
            j -= 1
        # If pointer is vertical, "consume" only a residue from the vertical sequence, and record a gap
        # in the other                                                                                 
        elif M[i][j].pointer == 'v':
            seq_1.append(a[i])
            seq_2.append('-')
            i -= 1
        # If pointer is horizontal, "consume" only a residue from the horizontal sequence, and record  
        # a gap in the other                                                                           
        elif M[i][j].pointer == 'h':
            seq_1.append('-')
            seq_2.append(b[j])
            j -= 1
    # Append final residues to sequences                                                               
    seq_1.append(a[i])
    seq_2.append(b[j])
    # Reverse sequence lists                                                                           
    seq_1.reverse()
    seq_2.reverse()
    # Transform sequence lists into strings                                                            
    str_1 = ''.join(seq_1)
    str_2 = ''.join(seq_2)
    # Return alignment strings                                                                         
    return str_1, str_2

# Needleman-Wunsch Global Alignment                                                                    
# Input: a and b, string sequences; subst, a substitution matrix; gap, negative gap penalty            
# Output: globally aligned sequence strings                                                            
def nw(a, b, subst, gap):
    # Generate scoring matrix                                                                          
    M = matrix(a, b, subst, gap)
    # Returned alignment strings                                                                       
    return backtrace(M, a, b)
