from Needleman_Wunsch import global_alignment
from helper_functions import parse

seq_1 = parse("seq3.txt")
seq_2 = parse("seq4.txt")

for alignment in global_alignment(seq_1, seq_2):
    print(alignment)
