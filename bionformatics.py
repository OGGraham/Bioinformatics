import helper_functions

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


class SmithWaterman:
    """
    1) Basic dynamic programming that runs in quadratic time and space [up to 50 marks].
     - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2):
        # Setup
        self.seq1, self.seq2, self.scoring_matrix, self.p = helper_functions.setup(seq1, seq2)

        # Penalty for matching with gap
        self.gap_penalty = -2

        # Setup both backtrack and scoring matrix
        self.scoring_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.scoring_matrix, local=True,
                                                                                   gap_penalty=self.gap_penalty)

    # Scoring function
    def score(self, a, b):
        # a & b match
        if a == b:
            return 3
        # Else
        return -3

    # Align 2 sequences
    def align(self):
        # Max score tracker
        max_score = -float('inf')
        max_index = []  # can have >1 greatest local alignment

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq[y-1], seq[x-1] as matrix has empty row & col at start
                    self.scoring_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.scoring_matrix[y-1][x] + self.gap_penalty,  # up
                    self.scoring_matrix[y][x - 1] + self.gap_penalty,  # left
                    0]  # 0 for local alignment
                # Update scoring matrix
                self.scoring_matrix[y][x] = max(vals)
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

        # Find all max alignments
        return helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)


class NeedlemanWunsch():
    """
    NeedlmanWunsch algorithm for global alignment (used in Hirschberg))
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2):

        # Setup
        self.seq1, self.seq2, self.scoring_matrix, self.p = helper_functions.setup(seq1, seq2)

        # Penalty for matching with gap
        self.gap_penalty = -1

        # Setup both backtrack and scoring matrix
        self.scoring_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.scoring_matrix, local=False,
                                                                                   gap_penalty=self.gap_penalty)

    # Scoring function
    def score(self, a, b):
        # a & b match
        if a == b:
            return 1
        # Else
        return -1

    # Align 2 sequences
    def align(self):
        # Global align -> always at bottom right index
        max_index = [len(self.backtrack_matrix)-1, len(self.backtrack_matrix[0])-1]

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq[y-1], seq[x-1] as matrix has empty row & col at start
                    self.scoring_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.scoring_matrix[y-1][x] + self.gap_penalty,  # up
                    self.scoring_matrix[y][x - 1] + self.gap_penalty,  # left
                ]
                # Update scoring matrix
                self.scoring_matrix[y][x] = max(vals)
                # Get index of max
                index = vals.index(max(vals))
                # Update backtrack matrix
                if index == 0:
                    self.backtrack_matrix[y][x] = 'D'
                elif index == 1:
                    self.backtrack_matrix[y][x] = 'U'
                elif index == 2:
                    self.backtrack_matrix[y][x] = 'L'


        # Find all max alignments
        return helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)


class Hirschberg():
    """
    2) Dynamic programming that runs in linear space [up to 65 marks].
     - for local alignment
     Implementation of Hirschberg's algorithm: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """
    def __init__(self, seq1, seq2):

        # Setup
        self.seq1, self.seq2, self.scoring_matrix, self.p = helper_functions.setup(seq1, seq2)

        # Penalty for matching with gap
        self.gap_penalty = -2

        # Setup both backtrack and scoring matrix
        self.scoring_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.scoring_matrix, local=True,
                                                                                   gap_penalty=self.gap_penalty)

        helper_functions.matrix_pretty_print(self.backtrack_matrix, self.seq1, self.seq2)

    # Substitution scoring function
    def substitute(self, a, b):
        # a & b match
        if a == b:
            return 2
        # Else
        return -1

    # Deletion scoring function
    def delete(self):
        return -2

    # Insert scoring function
    def insert(self):
        return -2

    # Last Row method (returns last row of scoring matrix - linear space complexity)
    def last_row(self, seq1, seq2):

        # print("Last Row Score on: {0} | {1}".format(seq1, seq2))

        # Init rows to 0s (as local alignment)
        prev_row = [0 for _ in range(len(seq2) + 1)]
        current_row = [0 for _ in range(len(seq2) + 1)]

        for j in range(1, len(seq2) + 1):
            prev_row[j] = prev_row[j-1] + self.insert()

        # Loop over seq2 and calc vals
        for i in range(1, len(seq1) + 1):
            current_row[0] = self.delete() + prev_row[0]
            for j in range(1, len(seq2) + 1):
                score_sub = prev_row[j - 1] + self.substitute(seq1[i - 1], seq2[j - 1])
                score_del = prev_row[j] + self.delete()
                score_ins = current_row[j - 1] + self.insert()
                current_row[j] = max(score_sub, score_del, score_ins)

            prev_row = current_row
            current_row = [0 for _ in range(len(seq2) + 1)]

            # print(prev_row)

        return prev_row

    # Calls recursive function with vals and returs
    def run(self):
        return self.align(self.seq1, self.seq2)

    # Hirschberg algorithm (ref. https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
    def align(self, seq1, seq2):
        out1, out2 = [], []

        # Empty seq1
        if len(seq1) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq2)):
                out1.append('-')
                out2.append(seq2[i])

        # Empty seq2
        elif len(seq2) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq1)):
                out2.append('-')
                out1.append(seq1[i])

        # Apply SW for optimal local alignment
        elif len(seq1) == 1 or len(seq2) == 1:
            NW = NeedlemanWunsch(seq1, seq2)
            out1, out2 = NW.align()

        else:
            # Get midpoint of Seq2
            seq2_mid = len(seq2) // 2

            # Get scoring of lhs (in linear space)
            r_left = self.last_row(seq2[:seq2_mid], seq1)
            # print(r_left)
            # Get scoring of rhs (in linear space) [reversed]
            r_right = self.last_row(seq2[seq2_mid:][::-1], seq1[::-1])
            r_right.reverse()
            # print(r_right)

            # Sum values and find argmax
            row = [l + r for l, r in zip(r_left, r_right)]
            maxidx, maxval = max(enumerate(row), key=lambda a: a[1])
            # print(row)

            # Partition seq1 at argmax
            seq1_mid = maxidx

            # Recursively call align on each half
            aligned_1_left, aligned_2_left = self.align(seq1[:seq1_mid], seq2[:seq2_mid])
            aligned_1_right, aligned_2_right = self.align(seq1[seq1_mid:], seq2[seq2_mid:])

            # Add results of recursive calls to  out
            out1 = aligned_1_left + aligned_1_right
            out2 = aligned_2_left + aligned_2_right


        print("IN:", seq1, seq2)
        print("OUT:", "".join(out1), "".join(out2))

        return out1, out2


class FASTA():
    """
    3) Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
    - for local alignment
    Implementation of the Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2):

        # Setup
        self.seq1, self.seq2, self.scoring_matrix, self.p = helper_functions.setup(seq1, seq2)

        # Setup both backtrack and scoring matrix
        self.scoring_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.scoring_matrix, local=True)

        # Penalty for matching with gap
        self.gap_penalty = -2

    def align(self):
        return [], []


if __name__ == "__main__":
    # Debug input 1 - example input from wiki (https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)
    # sequence1 = "TGTTACGG"  # seq1 = x
    # sequence2 = "GGTTGACTA"  # seq2 = y
    # Expected output: GG, TT, TT, -G, AA, CC

    # Debug input 2 - example from wiki (https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
    sequence1 = "TATGC"
    sequence2 = "AGTACGCA"

    # sequence1 = "TGTCC"
    # sequence2 = "ACTGACCT"


    print("Starting:")
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))
    print("------------")

    # Expected output: A-, G-, TT, AA, CT, GG, CC, A-
    # (w/ scoring del/ins = -2, sub(x,y) = +2 if match, else -1
    # exp1 = ['-', '-', 'T', 'A', 'T', 'G', 'C', '-']
    # exp2 = ['A', 'G', 'T', 'A', 'C', 'G', 'C', 'A']
    # print("Expected:")
    # print(helper_functions.alignment_pretty_print(exp1, exp2))

    # Part 1 - O(n^2) dynamic prog. (time + space)
    # results = part_one(p, scoring_matrix, seq1, seq2)
    # SW = SmithWaterman(sequence1, sequence2)
    # out1, out2 = SW.align()
    NW = NeedlemanWunsch(sequence1, sequence2)
    out1, out2 = NW.align()

    # Part 2 - O(n) dynamic prog. (space)
    # HB = Hirschberg(sequence1, sequence2)
    # out1, out2 = HB.run()

    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    # FA = FASTA
    # out1, out2 = FA.align()


    # Output - print results
    print("------------")
    print("Best Local Alignment:")
    helper_functions.alignment_pretty_print(out1, out2)

