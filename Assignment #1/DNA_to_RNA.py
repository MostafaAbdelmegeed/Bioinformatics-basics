from functions import parse, transcript

directory = input()

sequence = parse(directory)

print("DNA Sequence\t|\t{}\nRNA Sequence\t|\t{}\n".format(sequence, transcript(sequence)))
