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
     - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param p: string of all unique letters in seq1 & 2
    :param scoring_matrix: len(seq1) + 1 x len(seq2) + 1 matrix (+1 for the added - @ start of string)
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """
    # Init backtrack matrix
    backtrack_matrix = [[None for _ in range(len(scoring_matrix[0]))] for _ in range(len(scoring_matrix))]
    # Init first row & cols of matrices
    scoring_matrix[0] = [0 for _ in range(len(scoring_matrix[0]))]
    backtrack_matrix[0] = ['L' for _ in range(len(backtrack_matrix[0]))]
    for i in range(len(scoring_matrix)):
        scoring_matrix[i][0] = 0
        backtrack_matrix[i][0] = 'U'

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
    max_index = [[-1, -1]]  # can have >1 greatest local alignment

    # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there)
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
                max_index = [[y, x]]
            # Check if found another alignment with same max score
            elif max(vals) == max_score:
                max_index.append([y, x])

    # Find all max alignments
    alginments = []
    for item in max_index:
        # Start at max index & backtrack
        out = []
        x = item[1]
        y = item[0]
        while scoring_matrix[y][x] != 0:
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


if __name__ == "__main__":
    # Debug input - example input from wiki (https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)
    seq1 = "TGTTACGG"  # seq1 = x
    seq2 = "GGTTGACTA"  # seq2 = y
    p = str(set(seq1 + seq2))
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    scoring_matrix = [[None for x in range(len(seq1))] for y in range(len(seq2))]

    # Part 1)
    results = part_one(p, scoring_matrix, seq1, seq2)
    for item in results:
        alignment_pretty_print(item)
