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
            assert (len(a) == len(b))
            s = 0
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
        if not scoring_matrix or not indices:
            if len(a) == 3:
                scoring_matrix = self.scoring_matrix_3
                indices = self.scoring_matrix_3_indices
            elif len(a) == 2:
                scoring_matrix = self.scoring_matrix_2
                indices = self.scoring_matrix_2_indices
            else:
                scoring_matrix = self.scoring_matrix_1
                indices = self.scoring_matrix_1_indices
        return scoring_matrix[indices[a]][indices[b]]

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
    # Debug input 1
    alphabet = "ABC"
    scoring_matrix = [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
    sequence1 = "AABBAACA"
    sequence2 = "CBACCCBA"
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

    # Q2 - Part 2 - O(n^2) dynamic prog. (time + space)
    score, out1_indices, out2_indices, out1_chars, out2_chars = dynprogcost(sequence1, sequence2)

    # Output - print results
    print("------------")
    print("Score: {0}".format(score))
    print("Indices: {0} | {1}".format(out1_indices, out2_indices))
