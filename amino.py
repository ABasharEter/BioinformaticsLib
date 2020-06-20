from util import *
import re
from kmer import *
from collections import Counter
amino_map = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAU": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGU": "S",
    "AUA": "I",
    "AUC": "I",
    "AUG": "M",
    "AUU": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAU": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAU": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",
    "UAA": "",
    "UAC": "Y",
    "UAG": ".",
    "UAU": "Y",
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",
    "UGA": ".",
    "UGC": "C",
    "UGG": "W",
    "UGU": "C",
    "UUA": "L",
    "UUC": "F",
    "UUG": "L",
    "UUU": "F"
    }

iamino_map = {}
for k, v in amino_map.items():
    iamino_map.setdefault(v, []).append(k)

amino_mass = {
"G": 57,
"A": 71,
"S": 87,
"P": 97,
"V": 99,
"T": 101,
"C": 103,
"I": 113,
"L": 113,
"N": 114,
"D": 115,
"K": 128,
"Q": 128,
"E": 129,
"M": 131,
"H": 137,
"F": 147,
"R": 156,
"Y": 163,
"W": 186
}

def get_amino(rna):
    return "".join(amino_map["".join(mer)] for mer in chunks(rna,3))

def get_amino_mass(amino):
    return sum(amino_mass[a] if a.isalpha() else int(a) for a in amino)


def get_amino_rna(amino):
    if len(amino) == 1:
        return iamino_map[amino[0]]
    part2 = get_amino_rna(amino[1:])
    res = []
    for i in iamino_map[amino[0]]:
        for j in part2:
            print(i)
            print(j)
            res.append(i+j)
    return res

def find_peptide_encoding_dna2(dna, amino):
    rnas = set(get_amino_rna(amino))
    dnas = [rna_to_dna(rna) for rna in rnas] + [revers_complimernt(rna_to_dna(rna)) for rna in rnas]
    return [x for y in re.findall("|".join([f"(?=({x}))" for x in dnas]), dna) for x in y if len(x) > 0]
        
def search_once(dna, amino):
    res = []
    aminos = ""
    for i in range(0,len(dna)-2,3):
        mer = dna[i:i+3]
        rna = dna_to_rna(mer)
        aminos +=  amino_map[rna]
        if aminos[-len(amino):] == amino:
            res.append(dna[i-(len(amino)-1)*3:i+3])
    return res

def find_peptide_encoding_dna(dna, amino):
    res = search_once(dna, amino)
    res += search_once(dna[1:], amino)
    res += search_once(dna[2:], amino)
    res += [revers_complimernt(x) for x in search_once(revers_complimernt(dna), amino)]
    res += [revers_complimernt(x) for x in search_once(revers_complimernt(dna)[1:], amino)]
    res += [revers_complimernt(x) for x in search_once(revers_complimernt(dna)[2:], amino)]
    return res
    
def cyclospectrum(amino):
    res = [0]
    amino_double = amino+amino
    for l in range(1,len(amino)):
        for i in range(len(amino)):
            sub_amino = amino_double[i:i+l]
            res.append(get_amino_mass(sub_amino))
    res.append(get_amino_mass(amino))
    res = sorted(res)
    return res

def linearspectrum(amino):
    res = [0]
    for l in range(1,len(amino)):
        for i in range(len(amino)-l+1):
            sub_amino = amino[i:i+l]
            res.append(get_amino_mass(sub_amino))
    res.append(get_amino_mass(amino))
    res = sorted(res)
    return res

def amion_length_count(length):
    dp = [0]*(length+10)
    sizes = list(set(amino_mass.values()))
    for s in sizes:
        dp[s] += 1
    for i in range(length+10):
        dp[i] += sum(dp[i-s] for s in sizes if i-s >=0)
    return dp[length]

def cyclopeptide_sequencing(spectrom):
    candidates = set([""])
    final = []
    amion = set(amino_mass.keys())
    amion.remove("I")
    amion.remove("K")
    spectrom = sorted(spectrom)
    while len(candidates) > 0:
        candidates = set([c+a  for a in amion for c in candidates if all(pm in spectrom for pm in linearspectrum(c+a))])
        to_remove = []
        for c in candidates:
            if get_amino_mass(c) == spectrom[-1] and len(cyclospectrum(c)) == len(spectrom):
                final.append(c)
                to_remove.append(c)
        for r in to_remove:
            candidates.remove(r)
    return ["-".join([str(amino_mass[p]) for p in x]) for x in final]
                

def spectrom_score(spec, pep_spec):
    spec_list = list(spec)
    return sum((ps in spec_list,spec_list.remove(ps) if ps in spec_list else None)[0] for ps in pep_spec)

def leaderboard_cyclopeptide_sequencing(spectrom, n):
    candidates = set([("", 0)])
    leader_peptide = ""
    leader_score = 0
    amion = set(amino_mass.keys())
    amion.remove("I")
    amion.remove("K")
    spectrom = sorted(spectrom)
    while len(candidates) > 0:
        candidates = sorted([(c[0]+a,spectrom_score(spectrom,linearspectrum(c[0]+a)))  for a in amion for c in candidates],key=lambda  x: -x[1])
        to_remove = []
        if n < len(candidates):
            val = candidates[n-1][1]
            x = [i for i in range(n,len(candidates)) if candidates[i][1] >= val][-1]
            candidates = candidates[:x+1]
        for c,v in candidates:
            m = get_amino_mass(c)
            if m == spectrom[-1] and v > leader_score:
                leader_score = v
                leader_peptide = c
                to_remove.append((c,v))
            elif m > spectrom[-1]:
                to_remove.append((c,v))
        for r in to_remove:
            candidates.remove(r)
    return "-".join([str(amino_mass[p]) for p in leader_peptide])

def spectral_conv(spec):
    return [x-y for x in spec for y in spec if x-y >0]

def convolution_cyclopeptide_sequencing(spectrom, n, m):
    candidates = [([], 0)]
    leader_peptide = []
    leader_score = 0
    conv = spectral_conv(spectrom)
    amino = Counter([str(x) for x in conv if x >= 57 and x <= 200]).most_common()
    if len(amino) >= m:
        x = [i for i in range(m,len(amino)) if amino[i][1] >= amino[m][1]][-1]
        amino = amino[:x]
    amino = [x[0] for x in amino]
    amino += [str(get_amino_mass(x)) for x in amino_mass.keys()]
    amino = [[x] for x in set(amino)]
    spectrom = sorted(spectrom)
    while len(candidates) > 0:
        candidates = sorted([(c[0]+a,spectrom_score(spectrom,linearspectrum(c[0]+a))) for a in amino for c in candidates],key=lambda  x: -x[1])
        to_remove = []
        if n < len(candidates):
            val = candidates[n-1][1]
            x = [i for i in range(n,len(candidates)) if candidates[i][1] >= val][-1]
            candidates = candidates[:x+1]
        for c,v in candidates:
            m = get_amino_mass(c)
            if m == spectrom[-1] and v > leader_score:
                leader_score = v
                leader_peptide = c
                to_remove.append((c,v))
            elif m > spectrom[-1]:
                to_remove.append((c,v))
        for r in to_remove:
            candidates.remove(r)
    return "-".join(leader_peptide)
