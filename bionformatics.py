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

// TODO: For second part & introducing new scoring rules -> 'Gap Penalty (special penalty for consecutive “-”)'
//    https://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/03-05SequenceSimilarity.pdf (slide 31+)

"""

import helper_functions


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
    scoring_matrix, backtrack_matrix = helper_functions.matrix_setup(scoring_matrix)

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
    return helper_functions.backtrack(backtrack_matrix, max_indexes, seq1, seq2)


def part_two(p, scoring_matrix, seq1, seq2):
    """
    2) Dynamic programming that runs in linear space [up to 65 marks].
     - for local alignment
     Implementation of Hirschberg's algorithm: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    :param p: string of all unique letters in seq1 & 2
    :param scoring_matrix: len(seq1) + 1 x len(seq2) + 1 matrix (+1 for the added - @ start of string)
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """
    # Setup both backtrack and scoring matrix
    scoring_matrix, backtrack_matrix = helper_functions.matrix_setup(scoring_matrix)

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

    # NWScore method (returns last line of Needleman-Wunsch score matrix, Score(i, j))
    def NWScore(X, Y):

        # TODO: Score(0,0) = 0 // 2*length(Y) array - what does this mean???
        score_matrix = [[0 for _ in range(len(Y))] for _ in range(len(X))]
        raise NotImplementedError

        # TODO: See https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm for rest of implementation



    # Hirschberg algorithm (ref. https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
    def Hirschberg(X, Y):
        Z = ""
        W = ""
        # Empty X
        if len(X) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(1, len(Y)):
                Z += '-'
                W += Y[i]
        # Empty Y
        elif len(Y) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(X)):
                Z += X[i]
                W += '-'
        elif len(X) == 1 or len(Y) == 1:
            # TODO: NeedlemanWunsh -> think can just use Smith-Waterson from above
            # Z, W = NeedlemanWunsh(X, Y)
            raise NotImplementedError
        else:
            x_len = len(X)
            x_mid = x_len // 2
            y_len = len(Y)

            score_l = NWScore(X[:x_mid], Y)  # This should be a 1d arr
            score_r = NWScore(reversed(X[x_mid:]), reversed(Y))  # This should be a 1d arr
            # TODO: argmax of score_l + reversed(score_r) i.e. get max index
            # y_mid = arg_max(score_l + reversed(score_r))

            Z, W = Hirschberg(X[:x_mid], Y[:y_mid]) + Hirschberg(X[x_mid:], Y[y_mid:])

        return Z, W

    return []


if __name__ == "__main__":
    # Debug input - example input from wiki (https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)
    sequence1 = "AGTACGCA"  # seq1 = x
    sequence2 = "TATGC"  # seq2 = y

    # Setup
    seq1, seq2, scoring_matrix, p = helper_functions.setup(sequence1, sequence2)

    # Part 1 - O(n^2) dynamic prog. (time + space)
    # results = part_one(p, scoring_matrix, seq1, seq2)
    # Part 2 - O(n) dynamic prog. (space)
    results = part_two(p, scoring_matrix, seq1, seq2)
    # Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)

    # Output - print results
    print("Best Local Alignments:")
    for item in results:
        helper_functions.alignment_pretty_print(item)
