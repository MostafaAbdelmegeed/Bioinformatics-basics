# seq1.txt & seq2.txt are the sequences from slide #7 in Lecture_02
# The script will ask you for the sequences as pathes to the .txt files containing them
# simply just type seq1.txt and seq2.txt when asked for the sequences


import numpy as np



def parse(path):
    file = open(path, 'r')
    seq = file.read()
    seq = seq.rstrip('\r\n')
    return seq.upper()


H = 'H'
V = 'V'
D = 'D'
GAB = '_'


def inputs():
    print("Please enter First Sequence file path:")
    seq_1 = parse(input())
    print("Please enter Second Sequence file path:")
    seq_2 = parse(input())
    print('Please enter "Match" score: ')
    m_score = int(input())
    print('"Mis-match" score: ')
    mm_penalty = int(input())
    print('"Gap" score: ')
    g_penalty = int(input())
    return seq_1, seq_2, m_score, mm_penalty, g_penalty


def initialization(seq1, seq2, penalty=0):
    height, width = len(seq1), len(seq2)
    init_matrix = np.empty((height + 1, width + 1, 4), dtype=object)
    init_matrix[0, 0] = [0, None, GAB, GAB]
    for i in range(1, max(height + 1, width + 1)):
        if i < width + 1:
            init_matrix[0, i] = [i * penalty, H, GAB, seq2[i - 1]]
            init_matrix[1:, i, 3] = seq2[i - 1]
        if i < height + 1:
            init_matrix[i, 0] = [i * penalty, V, seq1[i - 1], GAB]
            init_matrix[i, 1:, 2] = seq1[i - 1]
    return init_matrix


def fill(matrix, mscore, mmpenalty, gpenalty):
    height, width = matrix.shape[0], matrix.shape[1]
    for i in range(1, height):
        for j in range(1, width):
            d = matrix[i - 1, j - 1, 0] + (mscore if matrix[i, j, 2] == matrix[i, j, 3] else mmpenalty)
            h = matrix[i, j - 1, 0] + gpenalty
            v = matrix[i - 1, j, 0] + gpenalty
            maximum = max(d, h, v)
            if maximum == d:
                if d == h:
                    matrix[i, j, 0:2] = [d, D + H]
                elif d == v:
                    matrix[i, j, 0:2] = [d, D + V]
                else:
                    matrix[i, j, 0:2] = [d, D]
            elif maximum == h:
                if h == d:
                    matrix[i, j, 0:2] = [h, D + H]
                if h == v:
                    matrix[i, j, 0:2] = [h, H + V]
                else:
                    matrix[i, j, 0:2] = [h, H]
            else:
                if v == d:
                    matrix[i, j, 0:2] = [v, D + V]
                elif v == h:
                    matrix[i, j, 0:2] = [v, H + V]
                else:
                    matrix[i, j, 0:2] = [v, V]
    return matrix


def traceback(matrix):
    height, width = matrix.shape[0], matrix.shape[1]
    sequence1 = ""
    sequence2 = ""
    index = (height - 1, width - 1)
    while index != (0, 0):
        if matrix[index[0], index[1], 1][0] == D:
            sequence1 += matrix[index[0], index[1], 2]
            sequence2 += matrix[index[0], index[1], 3]
            index = (index[0] - 1, index[1] - 1)
        elif matrix[index[0], index[1], 1][0] == H:
            sequence2 += matrix[index[0], index[1], 3]
            sequence1 += GAB
            index = (index[0], index[1] - 1)
        else:
            sequence1 += matrix[index[0], index[1], 2]
            sequence2 += GAB
            index = (index[0] - 1, index[1])
    return sequence1[::-1], sequence2[::-1]  # Inverting the strings


def global_alignment():
    seq1, seq2, match_score, mismatch_penalty, gap_penalty = inputs()
    matrix = initialization(seq1, seq2, penalty=gap_penalty)
    matrix = fill(matrix, match_score, mismatch_penalty, gap_penalty)
    sequence1, sequence2 = traceback(matrix)
    print("Best global alignment:\n{}\n{}".format(sequence1, sequence2))


global_alignment()
