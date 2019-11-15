import helper_functions
import re

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


class NeedlemanWunsch():
    """
    NeedlmanWunsch algorithm for global alignment (used in Hirschberg))
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Calculating score of pairings
        self.scoring_matrix = scoring_matrix

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(seq1, seq2)

        # Setup both backtrack and scoring matrix for global alignment
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=False,
                                                                                scoring_matrix=self.scoring_matrix,
                                                                                alphabet=self.alphabet,
                                                                                seq1=self.seq1,
                                                                                seq2=self.seq2)

    # Scoring function
    def score(self, a, b):
        # Get index in scoring matrix for chars a & b
        if a == '-':
            a_index = len(self.scoring_matrix[0]) - 1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix) - 1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Align 2 sequences
    def align(self):
        # Global align -> always at bottom right index
        max_index = [len(self.backtrack_matrix)-1, len(self.backtrack_matrix[0])-1]

        # Iterate over cost matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                vals = [
                    # seq[y-1], seq[x-1] as matrix has empty row & col at start
                    self.cost_matrix[y-1][x-1] + self.score(self.seq2[y-1], self.seq1[x-1]),  # diagonal
                    self.cost_matrix[y-1][x] + self.score(self.seq2[y-1], '-'),  # up
                    self.cost_matrix[y][x - 1] + self.score('-', self.seq1[x-1]),  # left
                ]
                # Update cost matrix
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

        # Score = value in bottom right cell (NW is global alignment)
        max_score = self.cost_matrix[len(self.cost_matrix)-1][len(self.cost_matrix[0])-1]

        # Find all max alignments
        return [max_score, helper_functions.backtrack(self.backtrack_matrix, max_index, self.seq1, self.seq2)]


class Hirschberg():
    """
    2) Dynamic programming that runs in linear space [up to 65 marks].
     - for local alignment
     Implementation of Hirschberg's algorithm: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2d array of each chars alignment
    """
    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Calculating score of pairings
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
            a_index = len(self.scoring_matrix[0]) - 1
        else:
            a_index = self.alphabet.index(a)
        if b == '-':
            b_index = len(self.scoring_matrix) - 1
        else:
            b_index = self.alphabet.index(b)

        # Return score from matrix
        return self.scoring_matrix[b_index][a_index]

    # Last Row method (returns last row of scoring matrix - linear space complexity)
    def last_row(self, seq1, seq2):
        max_val, max_index = -float('inf'), [1, 1]

        # Init rows to 0s (as local alignment)
        prev_row = [0 for _ in range(len(seq1) + 1)]
        current_row = [0 for _ in range(len(seq1) + 1)]

        # Init first row
        for j in range(1, len(seq1) + 1):
            prev_row[j] = max(0, prev_row[j-1] + self.score('-', seq1[j-1]))  # insert/left

        # Loop over remaining rows and calc scores
        for i in range(1, len(seq2) + 1):
            # Get first value in new row
            current_row[0] = max(0, prev_row[0] + self.score(seq2[i-1], '-'))  # del/up

            # Evaluate each value in row
            for j in range(1, len(seq1) + 1):
                score_sub = prev_row[j - 1] + self.score(seq2[i - 1], seq1[j - 1])  # diagonal
                score_ins = current_row[j - 1] + self.score('-', seq1[j-1])  # left
                score_del = prev_row[j] + self.score(seq2[i-1], '-')  # up

                # Local alignment -> max(vals, 0)
                current_row[j] = max(0, score_sub, score_del, score_ins)

                # Update max_val / index if score > max
                if current_row[j] > max_val:
                    max_val = current_row[j]
                    max_index = [i, j]  # y, x

            # Update prev row & clear current row
            prev_row = current_row
            current_row = [0 for _ in range(len(seq1) + 1)]

        return prev_row, max_val, max_index

    # Calls recursive function with vals and returs
    def run(self):
        # Get max index from forward pass
        _, max_val, max_index = self.last_row(self.seq1, self.seq2)

        # Get min index from backward pass
        _, _, min_index = self.last_row(self.seq1[::-1], self.seq2[::-1])

        # Subtract lengths from min index (s.t. actual start position)
        min_index[1] = len(self.seq1) - min_index[1]
        min_index[0] = len(self.seq2) - min_index[0]

        # Get local alignment
        return [self.align(self.seq1[min_index[1]:max_index[1]], self.seq2[min_index[0]:max_index[0]], min_index[1], min_index[0])]

    # TODO: need to track indices of sequences in parent string
    # Hirschberg algorithm (ref. https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
    def align(self, seq1, seq2, seq1_offset, seq2_offset):
        out1_chars, out2_chars = [], []
        out1_indices, out2_indices = [], []
        max_score = 0

        # Empty seq1
        if len(seq1) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq2)):
                out1_chars.append('-')
                out2_chars.append(seq2[i])
            # Score produced alignment
            prev_score = 0
            for i in range(len(out1_chars)):
                score = max(0, prev_score + self.score(out1_chars[i], out2_chars[i]))
                if score > max_score:
                    max_score = score
                prev_score = score

        # Empty seq2
        elif len(seq2) == 0:
            # Add -'s to remaining alignment s.t. valid
            for i in range(len(seq1)):
                out1_chars.append(seq1[i])
                out2_chars.append('-')
            # Score produced alignment
            prev_score = 0
            for i in range(len(out1_chars)):
                score = max(0, prev_score + self.score(out1_chars[i], out2_chars[i]))
                if score > max_score:
                    max_score = score
                prev_score = score

        # Apply SW for optimal local alignment
        elif len(seq1) == 1 or len(seq2) == 1:
            NW = NeedlemanWunsch(seq1, seq2, self.scoring_matrix, self.alphabet)
            needleman_output = NW.align()
            max_score, out1_indices, out2_indices, out1_chars, out2_chars = needleman_output[0], needleman_output[1][0], \
                                                                        needleman_output[1][1], needleman_output[1][2],\
                                                                        needleman_output[1][3]
            # TODO: adjust indices to account for offsets
            out1_indices = [x + seq1_offset for x in out1_indices]
            out2_indices = [x + seq2_offset for x in out2_indices]

        else:
            # Get midpoint of Seq2
            seq2_mid = len(seq2) // 2

            # Get scoring of lhs (in linear space)
            r_left, _, _ = self.last_row(seq1, seq2[:seq2_mid])

            # Get scoring of rhs (in linear space) [reversed]
            r_right, _, _ = self.last_row(seq1[::-1], seq2[seq2_mid:][::-1])
            r_right.reverse()  # flip back again for calculating total score

            # Sum values and find argmax
            row = [l + r for l, r in zip(r_left, r_right)]
            maxidx, maxval = max(enumerate(row), key=lambda a: a[1])

            # Partition seq1 at argmax
            seq1_mid = maxidx

            # Recursively call align on each half
            max_score_left, aligned_1_left_indices, aligned_2_left_indices, aligned_1_left_chars, aligned_2_left_chars = self.align(seq1[:seq1_mid], seq2[:seq2_mid], seq1_offset, seq2_offset)
            max_score_right, aligned_1_right_indices, aligned_2_right_indices, aligned_1_right_chars, aligned_2_right_chars = self.align(seq1[seq1_mid:], seq2[seq2_mid:], seq1_offset+seq1_mid, seq2_offset+seq2_mid)

            # Add results of recursive calls to out vars
            out1_chars = aligned_1_left_chars + aligned_1_right_chars
            out2_chars = aligned_2_left_chars + aligned_2_right_chars
            out1_indices = aligned_1_left_indices + aligned_1_right_indices
            out2_indices = aligned_2_left_indices + aligned_2_right_indices
            max_score = max_score_left + max_score_right

        # print("In: {0} | {1}".format(seq1, seq2))
        # print("Out: Score {0} - Indices {1} | {2} - Chars {3} | {4}".format(max_score, out1_indices, out2_indices, out1_chars, out2_chars))

        return max_score, out1_indices, out2_indices, out1_chars, out2_chars


# TODO: Use numba decorator to increase performance
class FASTA():
    """
    3) Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [up to 85 marks].
    - for local alignment
    :param seq1: sequence of chars, str
    :param seq2: sequence of chars, str
    :return: 2 arr's of each chars alignment
    """

    def __init__(self, seq1, seq2):
        # Setup
        self.seq1, self.seq2, self.scoring_matrix, self.p = helper_functions.create_cost_matrix(seq1, seq2)

        # Default word length
        self.word_length = 2

        # Penalty for matching with gap
        self.gap_penalty = -2

        # Setup both backtrack and scoring matrix
        self.scoring_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.scoring_matrix, local=True,
                                                                                   gap_penalty=self.gap_penalty)

    def seed(self, seq1, seq2):
        """
        Given 2 sequences and the word_length param, find the start and ends for all matching subwords in seq2 that
        are also in seq1
        :param seq1: str, seq1
        :param seq2: str, seq2
        :return: list of indexes of matching subwords
        """
        # TODO: find matches of len > word_length also?
        # Found words
        words = []
        # Get all words of length word_length and check if present in seq2
        for i in range(0, len(seq1)-self.word_length+1):
            # Get substring of word length
            word = seq1[i:i+self.word_length]
            # Get start & end indexes of matches
            matches = [(m.start(0), m.end(0)) for m in re.finditer(word, seq2)]
            if matches:
                for match in matches:
                    # Store in format (seq1 - start, end), (seq2 - start, end)
                    words.append([(i, i+self.word_length), match])
        return words

    def get_diagonals(self, words):
        """
        Given the word indices, find ones that lie on the same diagonal.
        They lie on the same diagonal if the difference in start index is the same
        :param words: list of tuples of  matches in words
        :return:
        """
        # Store the difference as the key, and the indices as the values
        diagonals = {}

        for item in words:
            # Get difference in starting index
            diff = item[0][0] - item[1][0]

            # Add item to dict
            try:
                diagonals[diff].append(item)
            except KeyError:
                diagonals[diff] = [item]

        return diagonals

    def combine(self, diagonals):
        """
        Given a dict of all seeds along the same diagonal, extend them if they lie upon the
        same diagonal and return the new start-end indices
        :param diagonals: dict, output of get_diagonals
        :return: 2d arr of tuples, [[(start1, end1), (start2, end2)], ...]
        """
        # Arr for output
        out = []

        # Combine the sequences that lie on the same diagonal
        for diagonal in diagonals.values():
            min_i1, max_i1 = float('inf'), 0  # seq1
            min_i2, max_i2 = float('inf'), 0  # seq2
            # Get min and max index of seeds on same diagonal
            for item in diagonal:
                if item[0][0] < min_i1:
                    min_i1 = item[0][0]
                if item[0][1] > max_i1:
                    max_i1 = item[0][1]
                if item[1][0] < min_i2:
                    min_i2 = item[1][0]
                if item[1][1] > max_i2:
                    max_i2 = item[1][1]

            # Produce new start-end indices for the sequences
            out.append([(min_i1, max_i1), (min_i2, max_i2)])

        return out

    def score(self, a, b):
        if a == b:
            return 2
        return -1

    def align(self):
        # 1) Seed - Find word sequences the sequences have in common, default size = 3
        # 2) Compile list of neighborhood words

        # Seed strings
        words = self.seed(self.seq1, self.seq2)
        # print(words)

        # Identify diagonals
        diagonals = self.get_diagonals(words)
        # print(diagonals)

        # Greedily extend words along diagonal as long as score improves
        words = self.combine(diagonals)
        print(words)

        # TODO: Run banded SmithWaterman on the greedily extended matches


        # # Iterate over words and score each one -> return the best alignment
        # max_score = -float('inf')
        # best_word = None
        #
        # for word in words:
        #     seq1_local = self.seq1[word[0][0], word[0][1]]
        #     seq2_local = self.seq2[word[1][0], word[1][1]]
        #     score = score(seq1_local, seq2_local)
        #     if score > max_score:
        #         max_score = score
        #         best_word = word
        #
        # # Return local alignment
        # return out1, out2


        exit(0)

        ## Max Score
        # max_score = 0
        # max1, max2 = [], []
        #
        # # Iterate over words and extend at either end
        # for word in words:
        #     # Get start & end indices
        #     start_seq1, end_seq1 = word[0][0], word[0][1]
        #     start_seq2, end_seq2 = word[1][0], word[1][1]
        #
        #     # Init score to matching score * word len
        #     score = self.score('a', 'a') * self.word_length
        #
        #     # For each match, extend in each direction until score = 0 (local)
        #     while score > 0 and (start_seq1 != 0 or start_seq2 != 0 or end_seq1 <= len(self.seq1)-1
        #                          or end_seq2 <= len(self.seq2)-1):
        #
        #         # Increment end pointers (if can)
        #         if start_seq1 > 0 and start_seq2 > 0:
        #             start_seq1 -= 1
        #             start_seq2 -= 1
        #             start_score = self.score(self.seq1[start_seq1], self.seq2[start_seq2])
        #             score = max(0, start_score + score)
        #             print("Start", start_seq1, start_seq2, start_score)
        #
        #         if end_seq1 <= len(self.seq1)-1 and end_seq2 <= len(self.seq2)-1:
        #             end_seq1 += 1
        #             end_seq2 += 1
        #             end_score = self.score(self.seq1[end_seq1], self.seq2[end_seq2])
        #             score = max(0, end_score + score)
        #             print("End", end_seq2, end_seq2, end_score)
        #
        #     # Store alignment if score > max
        #     if score > max_score:
        #         max_score = score
        #         max1, max2 = [x for x in self.seq1[start_seq1:end_seq1]],\
        #                      [x for x in self.seq2[start_seq2:end_seq2]]
        # return max1, max2


def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    SW = SmithWaterman(sequence1, sequence2, scoring_matrix, alphabet)
    results = SW.align()
    return results[0], results[1][0], results[1][1], results[1][2], results[1][3]


def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
    HB = Hirschberg(sequence1, sequence2, scoring_matrix, alphabet)
    results = HB.run()
    print("Output", results)
    return results[0][0], results[0][1], results[0][2], results[0][3], results[0][4]


if __name__ == "__main__":
    # # Debug input 1
    # alphabet = "ABC"
    # scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    # sequence1 = "AABBAACA"
    # sequence2 = "CBACCCBA"
    # Debug input 2
    alphabet = "ABCD"
    scoring_matrix = [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]]
    sequence1 = "AAAAACCDDCCDDAAAAACC"
    sequence2 = "CCAAADDAAAACCAAADDCCAAAA"
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
    # # Debug input 4
    # alphabet = "ABCD"
    # scoring_matrix = [
    #         [ 1,-5,-5,-5,-1],
    #         [-5, 1,-5,-5,-1],
    #         [-5,-5, 5,-5,-4],
    #         [-5,-5,-5, 6,-4],
    #         [-1,-1,-4,-4,-9]]
    # sequence1 = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
    # sequence2 = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
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
    # score, out1_indices, out2_indices, out1_chars, out2_chars = dynprog(alphabet, scoring_matrix, sequence1, sequence2)

    # Debug - O(n^2) dynamic prog. (time + space) -> [GLOBAL alignment]
    # NW = NeedlemanWunsch(sequence1, sequence2, scoring_matrix, alphabet)
    # results = NW.align()
    # score, out1_indices, out2_indices, out1_chars, out2_chars = results[0], results[1][0], results[1][1], results[1][2], results[1][3]

    # Part 2 - O(n) dynamic prog. (space)
    # FIXME: correct matches, incorrect scoring!
    score, out1_indices, out2_indices, out1_chars, out2_chars = dynproglin(alphabet, scoring_matrix, sequence1, sequence2)

    #  Part 3 - < O(n^2) heuristic procedure, similar to FASTA and BLAST (time)
    # FA = FASTA(sequence1, sequence2)
    # out1, out2 = FA.align()


    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
    # TODO: get this functional again
    # print("Best Local Alignment:")
    # helper_functions.alignment_pretty_print(out1_chars, out2_chars, sequence1, sequence2)


