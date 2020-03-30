from copy import deepcopy

import numpy as np

from functions import parse


def njm(dist_matrix):
    J_matrix = np.zeros(dist_matrix.shape)
    for i in range(J_matrix.shape[0]):
        for j in range(J_matrix.shape[1]):
            if j != i:
                J_matrix[i, j] = (
                        ((len(dist_matrix) - 2) * dist_matrix[i, j]) - np.sum(dist_matrix[i, :]) - np.sum(
                    dist_matrix[j, :]))
    return J_matrix


def njt(list_of_sequences):
    tree = []
    iteration = 0
    stack = []
    D_matrices = []
    J_matrices = []
    D_matrix = np.array([[0, 0.17, 0.59, 0.59, 0.77, 0.81, 0.87],
                         [0.17, 0, 0.6, 0.59, 0.77, 0.82, 0.86],
                         [0.59, 0.6, 0, 0.13, 0.75, 0.73, 0.86],
                         [0.59, 0.59, 0.13, 0, 0.75, 0.74, 0.88],
                         [0.77, 0.77, 0.75, 0.75, 0, 0.8, 0.93],
                         [0.81, 0.82, 0.73, 0.74, 0.8, 0, 0.9],
                         [0.87, 0.86, 0.86, 0.88, 0.93, 0.9, 0]])
    D_matrices.append(D_matrix.copy())
    while D_matrix.size > 4:
        J_matrix = njm(D_matrix).round(3)
        J_matrices.append(J_matrix)
        i, j = np.unravel_index(np.argmin(J_matrix, axis=None), J_matrix.shape)
        Dij = D_matrix[i, j]
        diu = (((1 / 2) * Dij) + ((1 / 2) * (1 / (len(D_matrix) - 2)) * (
                np.sum(D_matrix[i, :]) - np.sum(D_matrix[j, :])))).__round__(3)
        dju = (Dij - diu).__round__(3)
        stack.append([i, j, diu, dju])
        for k in range(len(D_matrix)):
            Dik = D_matrices[iteration][i, k]
            Djk = D_matrices[iteration][j, k]
            D_matrix[i, k] = ((1 / 2) * (Dik + Djk - Dij)).__round__(3)
        D_matrix = np.delete(D_matrix, j, axis=0)
        D_matrix = np.delete(D_matrix, j, axis=1)
        D_matrix[:, i] = D_matrix[i, :].T
        D_matrices.append(D_matrix.copy())
        if D_matrix.size == 4:
            indices = [0, 1, 2]
            if i > j:
                indices.pop(i)
                indices.pop(j)
            else:
                indices.pop(j)
                indices.pop(i)
            Dlast = (D_matrices[iteration][j, indices[0]] - diu).__round__(3)
            stack.append([1, Dlast])
        iteration += 1
    retrieve_absolute_indices(stack)
    return stack

## Needs to be implemented in recursively
def retrieve_absolute_indices(list_of_ids):
    absolute_indices = deepcopy(list_of_ids)
    shift = 0
    for i in range(1, len(absolute_indices) - 1):
        removed_index = absolute_indices[i - 1][1]
        if absolute_indices[i][0] >= removed_index:
            shift += 1
            for j in range(1, len(absolute_indices[1:])):
                absolute_indices[j][0] += shift
                absolute_indices[j][1] += shift
    print(list_of_ids)
    print(absolute_indices)


