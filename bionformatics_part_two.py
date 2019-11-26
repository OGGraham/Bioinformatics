import helper_functions
import re
import numpy as np
import math

"""
The second part is about coming up with crazy substitution scores. [total of 65 marks - the first question is for 
everyone and the second is for L4 only]
Design your own substitution-cost function that operates on pairs of sequences of letters instead of on pairs of 
letters. Clearly describe it on at most one page [15 marks]. 

For instance, such a function might give a cost of multi-letter substitution of "ABC" by "CBB", which is different 
from the simple addition of the single-letter costs, the mismatches AC and CB and the match BB;
have a fixed cost for inserting/deleting a sequence of Cs irrespective of its length, e.g. the first C incurs a 
cost of 1/2, the second one a cost of 1/4, the third one a cost of 1/8 and so on (so that the total is less than 1).

"""


class SmithWaterman2:
    """
    1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
     - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Scores of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=True)

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0])-1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix)-1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Align 2 sequences
    def align(self):
        # Max score tracker
        max_score = -float('inf')
        max_index = []  # can have >1 greatest local alignment

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                    self.cost_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.cost_matrix[y-1][x] + self.score(self.seq2[y-1], '-'),  # up
                    self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x-1]),  # left
                    0]  # 0 for local alignment
                # Update scoring matrix
                self.cost_matrix[y][x] = max(vals)
                # Get index of max
                index = vals.index(max(vals))
                # Update backtrack matrix
                if index == 0:
                    self.backtrack_matrix[y][x] = 'D'
                elif index == 1:
                    self.backtrack_matrix[y][x] = 'U'
                elif index == 2:
                    self.backtrack_matrix[y][x] = 'L'
                # Check if new greatest score seen
                if max(vals) > max_score:
                    max_score = max(vals)
                    max_index = [y, x]

        # Return max score and best local alignment chars + indices
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


def dynprogcost(sequence1, sequence2):
    # This section only expects alphabet of ABC
    if [x for x in set(sequence1 + sequence2) if x not in "ABC"]:
        print("Invalid sequence chars detected - expecting strings containing of ABC only.")
        exit(-1)

    SW = SmithWaterman2(sequence1, sequence2, scoring_matrix, alphabet)
    results = SW.align()
    return results[0], results[1][0], results[1][1], results[1][2], results[1][3]


if __name__ == "__main__":
    # # Debug input 1
    # alphabet = "ABC"
    # scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    # sequence1 = "AABBAACA"
    # sequence2 = "CBACCCBA"
    # # Debug input 2
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "AAAAACCDDCCDDAAAAACC"
    # sequence2 = "CCAAADDAAAACCAAADDCCAAAA"
    # # Debug input 3
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "AACAAADAAAACAADAADAAA"
    # sequence2 = "CDCDDD"
    # Debug input 4
    alphabet = "ABCD"
    scoring_matrix = [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]]
    sequence1 = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
    sequence2 = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # Q2 - Part 2 - O(n^2) dynamic prog. (time + space)
    score, out1_indices, out2_indices, out1_chars, out2_chars = dynprogcost(sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
