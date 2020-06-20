
dna_nucleotide_map = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

rna_nucleotide_map = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G"
}

dna_nucleotides = sorted(list(dna_nucleotide_map.keys()))
rna_nucleotides = sorted(list(rna_nucleotide_map.keys()))

def hamming(s1 ,s2):
    return sum(s1[i] != s2[i] for i in range(min(len(s1),len(s2)))) + abs(len(s1) - len(s2))

def revers_complimernt(gen, dna=True):
    return "".join(map(lambda x: dna_nucleotide_map[x] if dna else rna_nucleotide_map[x], gen[::-1]))

def dna_to_rna(dna):
    return dna.replace("T","U")
def rna_to_dna(dna):
    return dna.replace("U","T")


def enumerat_kmer(gen, k):
    for i in range(len(gen)-k+1):
        yield gen[i:k+i]

def generate_kmer(k, dna=True):
    if k == 0:
        yield ""
    else:
        for rest in generate_kmer(k-1, dna):
            for n in dna_nucleotides if dna else rna_nucleotides:
                yield n + rest

def gen_distance(gen, kmer):
    return sum(min(hamming(gen_kmer, kmer) for gen_kmer in enumerat_kmer(gen_str,len(kmer))) for gen_str in gen)
