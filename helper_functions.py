"""
Helper functions for rest of work
"""


def setup(seq1, seq2):
    # Generate set of all unique chars
    p = str(set(seq1 + seq2))
    # Generate scoring matrix
    scoring_matrix = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    # Return all
    return seq1, seq2, scoring_matrix, p


def matrix_setup(scoring_matrix, local, gap_penalty):
    """
    Given the scoring matrix, create the backtrack matrix & init the vals in both
    :param scoring_matrix: n x m array
    :param local: bool, whether local alignment or not
    :param gap_penalty: int (-ve), amount to penalize score by
    :return: n x m scoring matrix & n x m backtrack matrix (both w/ first row & col initialized)
    """
    # If local alignment, init scoring_matrix vals = 0; else = gap penalty
    if local:
        penalty = 0
    else:
        penalty = gap_penalty

    # Create backtrack matrix -> len + 1 as top left is blank
    backtrack_matrix = [[None for _ in range(len(scoring_matrix[0]))] for _ in range(len(scoring_matrix))]

    # Init first row & cols of matrices
    scoring_matrix[0] = [penalty * x for x in range(len(scoring_matrix[0]))]
    backtrack_matrix[0] = ['L' for _ in range(len(backtrack_matrix[0]))]
    for i in range(len(scoring_matrix)):
        scoring_matrix[i][0] = penalty * i
        backtrack_matrix[i][0] = 'U'

    # Set 0,0
    backtrack_matrix[0][0] = None

    return scoring_matrix, backtrack_matrix


def backtrack(backtrack_matrix, index, seq1, seq2):
    """
    Iterate over max_indexes and find all best local alignments
    :param backtrack_matrix: backtrack matrix
    :param index: arr of y, x coor of starting point for max alignment
    :param seq1: sequence str
    :param seq2: sequence str
    :return: 2d arr of all alignments
    """
    # Start at max index & backtrack
    out1, out2 = [], []
    # Stored as [y, x]
    x = index[1]
    y = index[0]
    # NB: seq1[x-1] & seq2[y-1] as backtrack & scoring matrix len = len(seq) + 1 (empty @ beginning)
    while backtrack_matrix[y][x] is not None and x >= 0 and y >= 0:
        if backtrack_matrix[y][x] == 'D':
            # Match both chars
            out1.append(seq1[x-1])
            out2.append(seq2[y-1])
            x -= 1
            y -= 1
        elif backtrack_matrix[y][x] == 'L':
            # Match seq1 w/ gap
            out1.append(seq1[x-1])
            out2.append('-')
            x -= 1
        else:
            # Match seq2 w/ gap
            out1.append('-')
            out2.append(seq2[y-1])
            y -= 1


    # Return alignment
    return list(reversed(out1)), list(reversed(out2))


def alignment_pretty_print(out1, out2):
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
    for i in range(len(out1)):
        s1.append(out1[i])
        s2.append(out2[i])
        if out1[i] == out2[i]:
            hyphens.append("|")
        else:
            hyphens.append(" ")
    # Print
    print(" ".join(s1))
    print(" ".join(hyphens))
    print(" ".join(s2))


def matrix_pretty_print(matrix, seq1, seq2):
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    print("   ", end="")
    for item in seq1:
        print(item, end="    ")
    print()

    for i in range(len(matrix)):
        print(seq2[i], matrix[i])
