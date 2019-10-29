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


def matrix_setup(scoring_matrix):
    """
    Given the scoring matrix, create the backtrack matrix & init the vals in both
    :param scoring_matrix: n x m array
    :return: n x m scoring matrix & n x m backtrack matrix (both w/ first row & col initialized)
    """
    # Create backtrack matrix
    backtrack_matrix = [[None for _ in range(len(scoring_matrix[0]))] for _ in range(len(scoring_matrix))]

    # Init first row & cols of matrices
    scoring_matrix[0] = [0 for _ in range(len(scoring_matrix[0]))]
    backtrack_matrix[0] = ['L' for _ in range(len(backtrack_matrix[0]))]
    for i in range(len(scoring_matrix)):
        scoring_matrix[i][0] = 0
        backtrack_matrix[i][0] = 'U'

    return scoring_matrix, backtrack_matrix


def backtrack(backtrack_matrix, max_indexes, seq1, seq2):
    """
    Iterate over max_indexes and find all best local alignments
    :param backtrack_matrix: backtrack matrix
    :param max_indexes: 2d arr of y, x coors of starting points for max alignments
    :param seq1: sequence str
    :param seq2: sequence str
    :return: 2d arr of all alignments
    """
    # Find all max alignments
    alginments = []
    for index in max_indexes:
        # Start at max index & backtrack
        out = []
        x = index[1]
        y = index[0]
        while backtrack_matrix[y][x] is not None:
            if backtrack_matrix[y][x] == 'D':
                # Match both chars
                out.append([seq1[x], seq2[y]])
                x -= 1
                y -= 1
            elif backtrack_matrix[y][x] == 'L':
                # Match seq1 w/ gap
                out.append([seq1[x], '-'])
                x -= 1
            else:
                # Match seq2 w/ gap
                out.append(['-', seq2[y]])
                y -= 1
        alginments.append(reversed(out))
    # Return all alignments
    return alginments


def part_one(p, scoring_matrix, seq1, seq2):
    """
    1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
     - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param p: string of all unique letters in seq1 & 2
    :param scoring_matrix: len(seq1) + 1 x len(seq2) + 1 matrix (+1 for the added - @ start of string)
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """
    # Setup both backtrack and scoring matrix
    scoring_matrix, backtrack_matrix = matrix_setup(scoring_matrix)

    # Scoring function
    def score(a, b):
        # a & b match
        if a == b:
            return 3
        # Else
        return -3
    # Penalty for matching with gap
    gap_penalty = -2

    # Max score tracker
    max_score = -float('inf')
    max_indexes = [[-1, -1]]  # can have >1 greatest local alignment

    # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
    for y in range(1, len(scoring_matrix)):  # y -> seq2
        for x in range(1, len(scoring_matrix[0])):  # x -> seq1
            vals = [
                scoring_matrix[y-1][x-1] + score(seq2[y], seq1[x]),  # diagonal
                scoring_matrix[y-1][x] + gap_penalty,  # up
                scoring_matrix[y][x - 1] + gap_penalty,  # left
                0]  # 0 for local alignment
            # Update scoring matrix
            scoring_matrix[y][x] = max(vals)
            # Get index of max
            index = vals.index(max(vals))
            # Update backtrack matrix
            if index == 0:
                backtrack_matrix[y][x] = 'D'
            elif index == 1:
                backtrack_matrix[y][x] = 'U'
            elif index == 2:
                backtrack_matrix[y][x] = 'L'
            # Check if new greatest score seen
            if max(vals) > max_score:
                max_score = max(vals)
                max_indexes = [[y, x]]
            # Check if found another alignment with same max score
            elif max(vals) == max_score:
                max_indexes.append([y, x])

    # Find all max alignments
    return backtrack(backtrack_matrix, max_indexes, seq1, seq2)


def alignment_pretty_print(alignment):
    """
    Given alignment, print in format:
        G C A T
        |   | |
        G T A T
    :param alignment: 2d arr of matches in alignment
    :return: Nothing
    """
    # Setup vars
    s1 = []
    hyphens = []
    s2 = []
    for item in alignment:
        s1.append(item[0])
        s2.append(item[1])
        if item[0] == item[1]:
            hyphens.append("|")
        else:
            hyphens.append(" ")
    # Print
    print(" ".join(s1))
    print(" ".join(hyphens))
    print(" ".join(s2))


def setup(seq1, seq2):
    # Generate set of all unique chars
    p = str(set(seq1 + seq2))
    # Add empty char onto front of strings
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    # Generate scoring matrix
    scoring_matrix = [[None for x in range(len(seq1))] for y in range(len(seq2))]
    # Return all
    return seq1, seq2, scoring_matrix, p


if __name__ == "__main__":
    # Debug input - example input from wiki (https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)
    sequence1 = "TGTTACGG"  # seq1 = x
    sequence2 = "GGTTGACTA"  # seq2 = y

    # Setup
    seq1, seq2, scoring_matrix, p = setup(sequence1, sequence2)

    # Part 1 - O(n^2) dynamic prog.
    results = part_one(p, scoring_matrix, seq1, seq2)

    # Output - print results
    print("Best Local Alignments:")
    for item in results:
        alignment_pretty_print(item)
