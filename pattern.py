import numpy as np 
import re

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
   
    return "".join(map(lambda x: nucleotide_map[x], gen[::-1]))

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

def motif_emnumeration(gen, k, d):
    p = set()
    for kmer in generate_kmer(k):
        if all(any(hamming(kmer, gen_kmer) <= d for gen_kmer in enumerat_kmer(gen_string, k)) for gen_string in gen):
            p.add(kmer)
    return list(p)

def gen_distance(gen, kmer):
    return sum(min(hamming(gen_kmer, kmer) for gen_kmer in enumerat_kmer(gen_str,len(kmer))) for gen_str in gen)


def median_string(gen, k):
    p = set()
    min_d = k+10
    for kmer in generate_kmer(k):
        d = gen_distance(gen, kmer)
        if d < min_d:
            min_d = d
            p = {kmer}
        elif d == min_d:
            p.add(kmer)
    return list(p)

def genom_matching(gen, pattern):
    return [m.start(1) for m in re.finditer(f"(?=({pattern}))", gen)]

def entropy(prop):
    return np.sum([-x*np.log2(x) for x in prop if x >0])

def motif_entropy(motifs, dna = True):
    n = len(motifs[0])
    result = 0
    for i in range(n):
        nucleotides = dna_nucleotides if dna else rna_nucleotides
        freq = {k:0 for k in nucleotides}
        for motif in motifs:
            freq[motif[i]] +=1
        tot = sum(list(freq.values()))
        result += entropy([v/tot for v in freq.values()])
    return result

def kmer_prob(kmer, prob):
    res = 1
    for i in range(len(kmer)):
        res *= prob[kmer[i]][i]
    return res

def most_prob_kmer(prob, gen, k):
    mx = -1
    mx_val = ""
    for kmer in enumerat_kmer(gen, k):
        p = kmer_prob(kmer, prob)
        if mx < p:
            mx = p
            mx_val = kmer
        elif mx == p:
            mx_val += " " + kmer 
    return mx_val


