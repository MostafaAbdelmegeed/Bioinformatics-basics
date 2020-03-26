import numpy as np

from distances import distances


def njm(list_of_sequences):
    dist_matrix = distances(list_of_sequences)
    J_matrix = np.zeros(dist_matrix.shape)

    for i in range(J_matrix.shape[0]):
        for j in range(J_matrix.shape[1]):
            if j != i:
                J_matrix[i, j] = (
                        ((len(list_of_sequences) - 2) * dist_matrix[i, j]) - np.sum(dist_matrix[i, :]) - np.sum(
                    dist_matrix[j, :]))
    return J_matrix
