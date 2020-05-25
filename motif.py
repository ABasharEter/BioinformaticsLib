import numpy as np 
import re
from kmer import *

def motif_emnumeration(gen, k, d):
    p = set()
    for kmer in generate_kmer(k):
        if all(any(hamming(kmer, gen_kmer) <= d for gen_kmer in enumerat_kmer(gen_string, k)) for gen_string in gen):
            p.add(kmer)
    return list(p)

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
    return mx_val


def motif_count(motifs, dna=True):
    nucleotides = dna_nucleotides if dna else rna_nucleotides
    counts = {n:[0]*len(motifs[0]) for n in nucleotides}
    for motif in motifs:
        for i in range(len(motif)):
            counts[motif[i]][i]+=1
    return counts

def motif_score(motifs, dna=True):
    counts = motif_count(motifs, dna)
    score = sum(
        len(motifs) - max([x[i] for x in counts.values()])
        for i in range(len(motifs[0]))
    )
    return score

def motif_profile(motifs, dna=True):
    counts = motif_count(motifs, dna)
    return {k:[v/len(motifs) for v in c] for k,c in counts.items()}
    
def motif_laplace_rule_count(motifs, dna=True):
    nucleotides = dna_nucleotides if dna else rna_nucleotides
    counts = {n:[1]*len(motifs[0]) for n in nucleotides}
    for motif in motifs:
        for i in range(len(motif)):
            counts[motif[i]][i]+=1
    return counts

def motif_laplace_rule_profile(motifs, dna=True):
    counts = motif_laplace_rule_count(motifs, dna)
    return {k:[v/(len(motifs) + 4) for v in c] for k,c in counts.items()}
        

def greedy_motif_search(gens, k, t):
    best_motifs = [g[:k] for g in gens]
    best_score = motif_score(best_motifs)
    for kmer in enumerat_kmer(gens[0], k):
        motifs = [kmer]
        for gen in gens[1:]:
            profile = motif_profile(motifs)
            next_kmer = most_prob_kmer(profile, gen, k)
            motifs.append(next_kmer)
        score = motif_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs

def greedy_motif_laplace_rule_search(gens, k, t):
    best_motifs = [g[:k] for g in gens]
    best_score = motif_score(best_motifs)
    for kmer in enumerat_kmer(gens[0], k):
        motifs = [kmer]
        for gen in gens[1:]:
            profile = motif_laplace_rule_profile(motifs)
            next_kmer = most_prob_kmer(profile, gen, k)
            motifs.append(next_kmer)
        score = motif_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs

def random_kmer(gen, k):
    s = np.random.randint(len(gen)-k)
    return gen[s:s+k]

def random_motif_search(gens, k, t):
    motifs = [random_kmer(gen,k) for gen in gens]
    best_motifs = motifs
    best_score = motif_score(best_motifs)
    while True:
        profile = motif_laplace_rule_profile(motifs)
        motifs = [most_prob_kmer(profile, gen, k) for gen in gens]
        score = motif_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
        else:
            return best_motifs



def profile_random_kmer(profile, gen, k):
    kmers = list(enumerat_kmer(gen, k))
    p = np.array([kmer_prob(kmer, profile) for kmer in kmers])
    p = p/np.sum(p)
    return np.random.choice(kmers, p = p)


def gibbs_samper(n, gens, k, t):
    motifs = [random_kmer(gen,k) for gen in gens]
    best_motifs = motifs
    best_score = motif_score(best_motifs)
    for _ in range(n):
        i = np.random.randint(t)
        profile = motif_laplace_rule_profile([motifs[j] for j in range(len(motifs)) if j != i])
        motifs[i] = profile_random_kmer(profile, gens[i], k)
        score = motif_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs

def run_n_times(func, times, **kwarg):
    best_motifs = func(**kwarg)
    best_score = motif_score(best_motifs)
    for _ in range(times):
        motifs = func(**kwarg)
        score = motif_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs