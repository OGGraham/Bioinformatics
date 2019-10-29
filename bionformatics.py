"""
The first part is about implementing different algorithms, exact and heuristic, for local alignment.
The input consists of four items: a string of all (unique) letters of length p, a (p + 1) x (p + 1) scoring matrix
(represented by a list of lists) and the two sequences.

You need to provide three functions which implement the following algorithms.

1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
2) Dynamic programming that runs in linear space [up to 65 marks].
3)A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].

They should output three items: the score of the best local alignment found by the algorithm plus two lists of indices,
one for each input sequences, that realise the matches/mismatches in the alignment.

This part will be marked automatically, so please make sure that you get both the input and the output in the right
format.

"""


def part_one(p, scoring_matrix, seq1, seq2):
    """
    1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
    Implementation of the Smith-Waterman algorithm
    :param p:
    :param scoring_matrix:
    :param seq1:
    :param seq2:
    :return:
    """
    # Init backtrack matrix
    backtrack_matrix = [[None for x in range(len(scoring_matrix[0]))] for y in range(len(scoring_matrix))]
    # Init first row & cols of matrices
    scoring_matrix[0] = [0 for x in range(len(scoring_matrix[0]))]
    backtrack_matrix[0] = ['L' for x in range(len(backtrack_matrix[0]))]
    for i in range(len(scoring_matrix)):
        scoring_matrix[i][0] = 0
        backtrack_matrix[i][0] = 'U'
    print(backtrack_matrix)

    # Define scoring function
    def score(a, b):
        # a & b match
        if a == b:
            return 3
        # Either a gap
        elif a != b and any(x == '-' for x in [a, b]):
            return -2
        # Just a mismatch
        else:
            return -3


    # Iterate over scoring matrix and generate scoring
    for i in range(1, len(scoring_matrix)):  # row
        for j in range(1, len(scoring_matrix[0])):  # col
            vals = [
                score(seq1[j - 1], seq2[i - 1]),  # diagonal
                score(seq1[j], seq2[i - 1]),  # left
                score(seq1[j - 1], seq2[i]),  # up
                0]
            # Update scoring matrix
            scoring_matrix[i][j] = max(vals)
            # Get index of max
            index = vals.index(max(vals))
            # Update backtrack matrix
            if index == 0:
                backtrack_matrix[i][j] = 'D'
            elif index == 1:
                backtrack_matrix[i][j] = 'L'
            elif index == 2:
                backtrack_matrix[i][j] = 'U'

    print(backtrack_matrix)





if __name__ == "__main__":
    # Debug input
    seq1 = "GTTAC"
    seq2 = "GTTGAC"
    p = str(set(seq1 + seq2))
    seq1 = " " + seq1
    seq2 = " " + seq2
    scoring_matrix = [[None for x in range(len(seq1))] for y in range(len(seq2))]

    # Part 1)
    part_one(p, scoring_matrix, seq1, seq2)
