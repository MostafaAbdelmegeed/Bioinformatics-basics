from functions import parse

print("Please enter Sequence file path:")
directory = input()

sequence = parse(directory)

c_count = sequence.count('C')
g_count = sequence.count('G')

print("Cytosine bases count\t|\t{0} Base\nGuanine bases count\t\t|\t{1} Base\nTotal number of bases\t|\t{2} Base\n"
      "CG content percentage\t|\t{3:3.2f}%\n".format(c_count, g_count, len(sequence),
                                                     (c_count + g_count) / len(sequence) * 100))
