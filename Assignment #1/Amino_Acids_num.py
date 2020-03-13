from functions import parse, translate

print("Please enter Sequence file path:")
directory = input()

sequence = parse(directory)

amino_acids, spared = translate(sequence)

print(
    "DNA Sequence\t\t\t|\t{}\nAmino Acids Present\t\t|\t{}\nAmino Acids count\t\t|\t{}\nSpared bases\t\t\t|\t{}\nCount of spared bases\t|\t{}\n".format(
        sequence,
        amino_acids,
        len(
            amino_acids),
        spared, len(spared)))
