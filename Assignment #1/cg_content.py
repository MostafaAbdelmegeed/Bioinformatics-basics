def parse(path):
    file = open(path, 'r')
    seq = file.read()
    seq = seq.rstrip('\r\n')
    return seq.upper()


sequence = parse("seq.txt")

c_count = sequence.count('C')
g_count = sequence.count('G')

print("Cytosine bases count\t|\t{0} Bases\nGuanine bases count\t\t|\t{1} Bases\nTotal number of bases\t|\t{2} Bases\n"
      "CG content percentage\t|\t{3:3.2f}%\n".format(c_count, g_count, len(sequence),
                                                     (c_count + g_count) / len(sequence) * 100))
