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
    :param p:
    :param scoring_matrix:
    :param seq1:
    :param seq2:
    :return:
    """
    pass


if __name__ == "__main__":
    # Debug input
    seq1 = "GTTAC"
    seq2 = "GTTGAC"
    p = str(set(seq1 + seq2))
    scoring_matrix = [[None for x in range(len(seq1)+1)] for y in range(len(seq2)+1)]

    # Part 1)
    part_one(p, scoring_matrix, seq1, seq2)
