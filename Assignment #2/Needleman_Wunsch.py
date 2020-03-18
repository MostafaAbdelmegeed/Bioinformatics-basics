# seq1.txt & seq2.txt are the sequences from slide #7 in Lecture_02
# The script will ask you for the sequences as paths to the .txt files containing them
# simply just type seq1.txt and seq2.txt when asked for the sequences


import numpy as np


def parse(path):  # Parse .txt files for sequence
    file = open(path, 'r')
    seq = file.read()
    seq = seq.rstrip('\r\n')
    return seq.upper()


D = 'D'  # Diagonal Direction
V = 'V'  # Vertical Direction
H = 'H'  # Horizontal Direction
DIRECTIONS = [V, H, D]
GAB = '_'
PERMUTATIONS = []  # Indices of slots with permutation potential


def figure_direction(maximum, values):  # Returns all the directions possible for the currents slot
    direction = ''
    if all(value == maximum for value in values):
        direction = V + H + D
    else:
        for i in range(len(values)):
            if values[i] == maximum:
                direction += DIRECTIONS[i]
    return direction


def get_permutations(matrix):  # Returns all possible Alignments with maximum score
    perm_matrix = np.copy(matrix)
    all_permutations = []
    while PERMUTATIONS:
        perm = PERMUTATIONS.pop(0)
        perm_matrix[perm[0], perm[1], 1] = perm[2]
        all_permutations.append(traceback(perm_matrix))
    return all_permutations


# Returns a 3D matrix with first row and first column initial values and 2 matrices containing the protein bases
# [:,:,0] Contains the values of the matrix
# [:,:,1] Contains the directions of the matrix
# [:,:,2] Contains the first sequence protein bases --- as a matrix of repeated rows
# [:,:,3] Contains the second sequence protein bases --- as a matrix of repeated columns
def initialization(seq1, seq2, penalty=0):
    width, height = len(seq1), len(seq2)
    init_matrix = np.empty((height + 1, width + 1, 4), dtype=object)
    init_matrix[0, 0] = [0, None, GAB, GAB]
    for i in range(1, max(height + 1, width + 1)):
        if i < width + 1:
            init_matrix[0, i] = [i * penalty, H, GAB, seq1[i - 1]]
            init_matrix[1:, i, 3] = seq1[i - 1]
        if i < height + 1:
            init_matrix[i, 0] = [i * penalty, V, seq2[i - 1], GAB]
            init_matrix[i, 1:, 2] = seq2[i - 1]
    return init_matrix


def fill(matrix, mscore, mmpenalty, gpenalty):  # Returns a filled matrix with values according to the score given
    height, width = matrix.shape[0], matrix.shape[1]
    for i in range(1, height):
        for j in range(1, width):
            d = matrix[i - 1, j - 1, 0] + (mscore if matrix[i, j, 2] == matrix[i, j, 3] else mmpenalty)
            v = matrix[i - 1, j, 0] + gpenalty
            h = matrix[i, j - 1, 0] + gpenalty
            maximum = max(v, h, d)
            directions = figure_direction(maximum, [v, h, d])
            matrix[i, j, 0:2] = [maximum, directions]
    return matrix


def traceback(matrix):
    height, width = matrix.shape[0], matrix.shape[1]
    sequence1 = ""
    sequence2 = ""
    index = (height - 1, width - 1)
    while index != (0, 0):
        directions = matrix[index[0], index[1], 1]  # Directions possible for the current slot
        if len(directions) > 1:  # Possible permutations
            for each_direction in directions[
                                  1:]:  # Add all the directions except the first as it's already getting it's own run through this iteration
                PERMUTATIONS.append((index[0], index[1], each_direction))
        if directions[0] == D:  # If direction is Diagonal
            sequence2 += matrix[index[0], index[1], 2]
            sequence1 += matrix[index[0], index[1], 3]
            index = (index[0] - 1, index[1] - 1)  # Move on the Diagonal
        elif directions[0] == H:  # If direction is Horizontal
            sequence1 += matrix[index[0], index[1], 3]
            sequence2 += GAB  # Base deletion
            index = (index[0], index[1] - 1)  # Move on the Horizontal
        else:  # Else Vertical
            sequence2 += matrix[index[0], index[1], 2]
            sequence1 += GAB  # Base deletion
            index = (index[0] - 1, index[1])  # Move on the Vertical
    return sequence1[::-1], sequence2[::-1]  # Inverting the strings


def global_alignment(seq1, seq2, match_score, mismatch_penalty, gap_penalty, perm=False):
    matrix = initialization(seq1, seq2, penalty=gap_penalty)  # Check Function Header
    matrix = fill(matrix, match_score, mismatch_penalty, gap_penalty)  # Check Function Header
    if not perm:  # No Permutations, just the first alignment
        alignment = traceback(matrix)
        print("Possible Alignment with score {}\n{}\n{}".format(matrix[-1, -1, 0], alignment[0], alignment[1]))
    else:  # All Possible Alignments
        alignments = [(traceback(matrix))]
        permutations = get_permutations(matrix)
        for perm in permutations:
            alignments.append(perm)
        for alignment in alignments:
            print("Possible Alignment with score {}\n{}\n{}".format(matrix[-1, -1, 0], alignment[0], alignment[1]))
            print('---------------')


print("Please Enter First Sequence Path: ex. 'seq1.txt':")
seq1_path = input()
print("Please Enter Second Sequence Path: ex. 'seq2.txt':")
seq2_path = input()
print("Please Enter Matching Score:")
match_score = int(input())
print("Please Enter Mismatch Penalty:")
mismatch_penalty = int(input())
print("Please Enter Gap Penalty:")
gap_penalty = int(input())
print("Do you want all the permutations? (Y/N)\tThis may take some time with long sequences")
perm = True if input().upper() == 'Y' else False
global_alignment(parse(seq1_path), parse(seq2_path), match_score, mismatch_penalty, gap_penalty, perm=perm)
