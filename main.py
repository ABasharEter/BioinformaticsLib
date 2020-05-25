from motif import *

if __name__ == "__main__":
    gens = """TATGAGGTC
GCCCTAGA
AAATAGAT
TTGTGCTA""".split("\n")
    motifs = """GTC
CCC
ATA
GCT""".split("\n")
    best_motifs = motifs
    profile = motif_laplace_rule_profile(motifs)
    motifs = [most_prob_kmer(profile, gen, 3) for gen in gens]
    print(" ".join(motifs))
