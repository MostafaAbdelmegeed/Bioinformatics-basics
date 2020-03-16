THYMINE = 'T'
URACIL = 'U'


def parse(path):
    file = open(path, 'r')
    seq = file.read()
    seq = seq.rstrip('\r\n')
    return seq.upper()


def transcript(template_strand):
    return template_strand.replace(THYMINE, URACIL)


sequence = parse("seq.txt")

print("DNA Sequence\t|\t{}\nRNA Sequence\t|\t{}\n".format(sequence, transcript(sequence)))
