import numpy as np

from Needleman_Wunsch import global_alignment


def distances(list_of_sequences, used_align=False):
    distances_matrix = np.zeros((len(list_of_sequences), len(list_of_sequences)))
    used_alignments = []
    for i in range(len(list_of_sequences)):
        for j in range(len(list_of_sequences)):
            alignment = global_alignment(list_of_sequences[i], list_of_sequences[j])
            mismatch_count = 0
            gap_count = 0
            for k in range(len(alignment[0])):
                if alignment[0][k] == '_' or alignment[1][k] == '_':
                    gap_count += 1
                elif alignment[0][k] != alignment[1][k]:
                    mismatch_count += 1
            used_alignments.append(alignment)
            distances_matrix[i, j] = gap_count / (len(alignment[0]) - gap_count)
    if used_align:
        return distances_matrix, used_alignments
    else:
        return distances_matrix
