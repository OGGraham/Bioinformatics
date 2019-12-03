import helper_functions
import itertools


"""
The second part is about coming up with crazy substitution scores. [total of 65 marks - the first question is for 
everyone and the second is for L4 only]
Design your own substitution-cost function that operates on pairs of sequences of letters instead of on pairs of 
letters. Clearly describe it on at most one page [15 marks]. 

Stefan's alternative:

It is to solve 2 above with this specific function: it is defined by a symmetric matrix which is indexed by all single 
letters plus some short substrings. For instance, it could have ABC and CBB as indices and the cost of the "mismatch" 
could be different from the sum of the scores of the mismatches AC and CB and the match BB.

Scoring function:

- Generate scores for subsequences:
    - If they dont contain any -'s AND >= 50% of chars match -> score = sum of individual char matches + len subseq
    - Else, 0
    
- For scoring sequences/backtracking:
    - Iterate over sequence in groups of chars len 3
    - E.g. AAABC -> AAA, AAB, ABC

"""


class SmithWaterman2:
    def __init__(self, seq1, seq2, scoring_matrix, alphabet):
        # Part 2 - Setup multiple scoring matrices for subsequences
        self.scoring_matrix_1, self.scoring_matrix_1_indices,\
        self.scoring_matrix_2, self.scoring_matrix_2_indices,\
        self.scoring_matrix_3, self.scoring_matrix_3_indices, = self.scoring_matricies_setup(scoring_matrix, alphabet)

        # Set of unique characters (same order as in scoring matrix)
        self.alphabet = alphabet

        # Sequences
        self.seq1 = seq1
        self.seq2 = seq2

        # Setup cost matrix
        self.cost_matrix = helper_functions.create_cost_matrix(self.seq1, self.seq2)

        # Setup both backtrack and cost matrix initial values
        self.cost_matrix, self.backtrack_matrix = helper_functions.matrix_setup(self.cost_matrix, local=True)

    # Setup cost matrices for subsequence scoring
    def scoring_matricies_setup(self, scoring_matrix, alphabet):

        # Local function for calculating scores of subsequences
        def get_scores(a, b, indicies):
            # Count matches in a vs. b
            matches = 0
            for i in range(len(a)):
                if a[i] == b[i]:
                    matches += 1

            # Only calc new score for seqs not containing -'s and >= 50% of matching chars
            if "-" in set(a+b) or matches < len(a) * 0.5:
                return None
            else:
                # Sequences get a scoring bonus of their length
                s = len(a)
                # Score += sum of matching each char against each other
                for i in range(len(a)):
                    s += self.score(a[i], b[i], scoring_matrix, indicies)
                return s

        # Local function for creating scoring matrix
        def create_scoring_matrix(new_indices, indicies):
            # Create n x n size matrix
            matrix = [[0 for _ in range(len(new_indices))] for _ in range(len(new_indices))]
            # Populate w/ correct scores
            for x, x_val in new_indices.items():
                for y, y_val in new_indices.items():
                    matrix[y_val][x_val] = get_scores(x, y, indicies)
            return matrix

        # Get scoring matrix & indices for len 1
        matrix_len_1 = scoring_matrix
        matrix_len_1_indices = {v: k for k, v in enumerate(alphabet + '-')}
        # Get scoring matrix & indices for len 2
        matrix_len_2_indices = {v: k for k, v in enumerate([''.join(x) for x in itertools.product(alphabet + '-',
                                                                                                   repeat=2)])}
        matrix_len_2 = create_scoring_matrix(matrix_len_2_indices, matrix_len_1_indices)
        # Get scoring matrix & indices for len 3
        matrix_len_3_indices = {v: k for k, v in enumerate([''.join(x) for x in itertools.product(alphabet + '-',
                                                                                                   repeat=3)])}
        matrix_len_3 = create_scoring_matrix(matrix_len_3_indices, matrix_len_1_indices)
        return matrix_len_1, matrix_len_1_indices, matrix_len_2, matrix_len_2_indices, matrix_len_3, matrix_len_3_indices

    # Scoring function
    def score(self, a, b, scoring_matrix=None, indices=None):

        # Local functions for scoring a, b of different lengths:
        def score_len_3(a, b):
            # Check if score for len(3) subseqs
            id1, id2 = self.scoring_matrix_3_indices[a], self.scoring_matrix_3_indices[b]
            s = self.scoring_matrix_3[id1][id2]
            if s:
                return s
            # Get best score for len(2) subseqs + the other match
            s = []
            a_subseqs = [[a[:2], a[2]], [a[1:], a[0]]]
            b_subseqs = [[b[:2], b[2]], [b[1:], b[0]]]
            for a_item in a_subseqs:
                id1 = self.scoring_matrix_2_indices[a_item[0]]
                for b_item in b_subseqs:
                    id2 = self.scoring_matrix_2_indices[b_item[0]]
                    v = self.scoring_matrix_2[id1][id2]
                    if v is not None:
                        v += self.scoring_matrix_1[self.scoring_matrix_1_indices[a_item[1]]][
                            self.scoring_matrix_1_indices[b_item[1]]]
                        s.append(v)
            if any(s):
                return max(s)
            # Get score of individual matching chars
            s = 0
            for i in range(len(a)):
                s += self.scoring_matrix_1[self.scoring_matrix_1_indices[a[i]]][self.scoring_matrix_1_indices[b[i]]]
            return s

        def score_len_2(a, b):
            s = self.scoring_matrix_2[self.scoring_matrix_2_indices[a]][self.scoring_matrix_2_indices[b]]
            if s is not None:
                return s
            else:
                s = 0
                for i in range(len(a)):
                    s += self.scoring_matrix_1[self.scoring_matrix_1_indices[a[i]]][self.scoring_matrix_1_indices[b[i]]]
                return s

        def score_len_1(a, b):
            return self.scoring_matrix_1[self.scoring_matrix_1_indices[a]][self.scoring_matrix_1_indices[b]]

        # Input just a & b -> 1 <= len(a & b) <= 3
        if not scoring_matrix or not indices:
            assert (len(a) == len(b))
            if len(a) == 3:
                return score_len_3(a, b)
            elif len(a) == 2:
                return score_len_2(a, b)
            else:
                return score_len_1(a, b)
        # Input a, b, scoring matrix & indices
        else:
            return scoring_matrix[indices[a]][indices[b]]

    # Partial backtrack to obtain subsequences for scoring
    def get_subsequences(self, x, y):
        """
        Given starting coordinate of x, y (square finding score for), return the subsequences of len <= 3 appropriately
        from that starting position
        :param x: x coor, int
        :param y: y coor, int
        :return: 3 appropriate susequences for scoring -> diagonal, up, left
        """
        # Local function - partial backtrack (maximum of 2 extra chars)
        def partial_backtrack(x_current, y_current):
            current_score = self.cost_matrix[y_current][x_current]
            subseq1, subseq2 = "", ""
            while current_score != 0 and len(subseq1) < 2:
                # Add chars into subseq
                subseq1 += self.seq1[x_current - 1]
                subseq2 += self.seq2[y_current - 1]
                # Move in direction of best local alignment
                score_diag = self.cost_matrix[y_current - 1][x_current - 1] if y_current-1 >= 0 and x_current-1 >= 0 else 0
                score_left = self.cost_matrix[y_current][x_current - 1] if x_current-1 >= 0 else 0
                score_up = self.cost_matrix[y_current - 1][x_current] if y_current-1 >= 0 else 0
                max_v = max(score_diag, score_left, score_up)
                # Move in dir of best score
                if max_v == score_diag:
                    x_current -= 1
                    y_current -= 1
                elif max_v == score_up:
                    y_current -= 1
                else:
                    x_current -= 1
                # Update current score (of cell looking at)
                current_score = self.cost_matrix[y_current][x_current]
            return subseq1, subseq2

        # Init subseqs to their starting chars
        diag1, diag2 = self.seq1[x-1], self.seq2[y-1]
        up1, up2 = '-', self.seq2[y-1]
        left1, left2 = self.seq1[x-1], '-'

        # Diagonal
        r1, r2 = partial_backtrack(x-1, y-1)
        diag1 += r1
        diag2 += r2
        # Up
        r1, r2 = partial_backtrack(x, y-1)
        up1 += r1
        up2 += r2
        # Left
        r1, r2 = partial_backtrack(x-1, y)
        left1 += r1
        left2 += r2

        return [[diag1, diag2], [up1, up2], [left1, left2]]

    # Align 2 sequences
    def align(self):
        # Max score tracker
        max_score = -float('inf')
        max_index = []  # can have >1 greatest local alignment

        # Iterate over scoring matrix and generate scoring (start at 1,1 and work from there) (O(n^2) method)
        for y in range(1, len(self.seq2)+1):  # y -> seq2
            for x in range(1, len(self.seq1)+1):  # x -> seq1
                # Get subsequences for scoring (partial backtrack)
                subseqs = self.get_subsequences(x, y)
                vals = [
                    # seq2[y-1], seq1[x-1] as matrix has empty row & col at start
                    self.score(subseqs[0][0], subseqs[0][1]) + self.cost_matrix[y-1][x-1],  # diagonal
                    self.score(subseqs[1][0], subseqs[1][1]) + self.cost_matrix[y-1][x],  # up
                    self.score(subseqs[2][0], subseqs[2][1]) + self.cost_matrix[y][x-1],  # left
                    0]  # 0 for local alignment
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
    if [x for x in set(sequence1 + sequence2) if x not in {"A", "B", "C"}]:
        print("Invalid sequence chars detected - expecting strings containing ABC only.")
        exit(-1)
    # Default params -> scoring matrix and ABC alphabet
    scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    alphabet = "ABC"
    SW = SmithWaterman2(sequence1, sequence2, scoring_matrix, alphabet)
    results = SW.align()
    return results[0], results[1][0], results[1][1]


if __name__ == "__main__":
    # # Debug input 1
    # sequence1 = "AABBAACA"
    # sequence2 = "CBACCCBA"
    # Test input 2
    sequence1 = "AAAABCABABCAABCBA"
    sequence2 = "BBABAAABCCCBABCAA"

    print("Starting:")
    # Strip to ensure no whitespace
    sequence1, sequence2 = sequence1.strip(), sequence2.strip()
    print("Seq 1 - {0} ".format(sequence1))
    print("Seq 2 - {0}".format(sequence2))

    # Q2 - Part 2 - O(n^2) dynamic prog. (time + space)
    score, out1_indices, out2_indices = dynprogcost(sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
