"""
Helper functions for rest of work
"""


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
