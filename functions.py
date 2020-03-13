ADENINE = 'A'
CYTOSINE = 'C'
GUANINE = 'G'
THYMINE = 'T'
URACIL = 'U'

CODONS = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
    'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
}


def parse(path):
    file = open(path, 'r')
    seq = file.read()
    seq = seq.rstrip('\r\n')
    return seq


def transcript(template_strand):
    return template_strand.replace(THYMINE, URACIL)


def translate(seq):
    transcripted_strand = transcript(seq)
    transcripted_strand, spare = insure_divisability(transcripted_strand, 3)
    polypeptide = ""
    for i in range(0, len(transcripted_strand), 3):
        polypeptide += CODONS[transcripted_strand[i:i + 3]]
    return polypeptide, spare


def insure_divisability(seq, factor=3):
    spare = len(seq) % factor
    if spare == 0:
        return seq, None
    elif spare == 1:
        base = seq[len(seq) - 1]
        seq = seq[:-1]
        return seq, base
    elif spare == 2:
        bases = seq[len(seq) - 2:]
        seq = seq[:-2]
        return seq, bases


def amino_count(seq):
    transcription = transcript(seq)
    polypeptide, spare = translate(transcription)
    return len(polypeptide)
