from motif import *
from collections import defaultdict


def composition_k(gen, k):
    return sorted(list(enumerat_kmer(gen, k)))

def path_to_genome(path):
    return path[0][:-1] + "".join((x[-1] for x in path))

def overlapping_graph(kmers):
    perfixs = defaultdict(list)
    for kmer in kmers:
        perfixs[kmer[:-1]].append(kmer)
    suffixs = defaultdict(list)
    for kmer in kmers:
        for per in perfixs[kmer[1:]]:
            suffixs[kmer].append(per)
    return suffixs

def de_bruijn_graph(gen, k):
    k = k-1
    res = defaultdict(list)
    for i in range(len(gen) - k):
        res[gen[i:i+k]].append(gen[i+1:i+k+1])
    return res

def de_bruijn_graph_kmers(kmers):
    suffixs = defaultdict(set)
    for kmer in kmers:
        suffixs[kmer].add(kmer[1:])
    perfixs = defaultdict(list)
    for kmer in kmers:
        for suf in suffixs[kmer]:
            perfixs[kmer[:-1]].append(suf)
            if suf not in perfixs:
                perfixs[suf] = []
    return perfixs

def paired_de_bruijn_graph_kmers(kmers):
    g = defaultdict(list)
    for kmer in kmers:
        e1,e2 = kmer.split("|")
        v1 = f"{e1[:-1]}|{e2[:-1]}"
        v2 = f"{e1[1:]}|{e2[1:]}"
        g[v1].append(v2)
        if v2 not in g:
            g[v2] = []
    return g

def constrac_paired_string(path, d ,k):
    x = path_to_genome(list(x.split("|")[0] for x in path))
    y = path_to_genome(list(x.split("|")[1] for x in path))
    return x + y[-d-k:]

def solve_paired_string(g, d ,k):
    e_count = sum(len(x) for x in g.values())
    print(e_count)
    for n in g.keys():
        path = eularian(g, n)
        x = path_to_genome(list(x.split("|")[0] for x in path))
        y = path_to_genome(list(x.split("|")[1] for x in path))
        p1 = y[:e_count-1]
        p2 = x[-(e_count-1):]
        print(x)
        print(y)
        print(p1==p2)
        sol = x + y[-d-k:]
    if len(sol) == e_count+k-1+k+d:
        return sol
    return ""

def inverse_graph(g):
    res = defaultdict(list)
    for v,adj in g.items():
        for u in adj:
            res[u].append(v)
        if v not in res:
            res[v] = []
    return res

def eularian(g, start=None):
    stack = []
    ig = inverse_graph(g)
    if start is not None:
        stack.append(start)
    else:
        v = [x for x in g.keys() if len(g[x]) > len(ig[x])]
        if v:
            stack.append(v[0])
        else:
            stack.append(list(g.keys())[0])
    res = []
    adj = {k:list(v) for k,v in g.items()}
    while len(stack)>0:
        v = stack[-1]
        if len(adj[v]) == 0:
            res.append(v)
            stack.pop()
        else:
            u = adj[v].pop()
            stack.append(u)
    return res[::-1]

def parse_graph(data):
    res = {}
    for line in data.split("\n"):
        k,v = line.split("->")
        v = v.split(",")
        k = k.strip()
        res[k] = [x.strip() for x in v]
        for x in v:
            x = x.strip()
            if x not in res:
                res[x] = []
    return res


def gen_universal_strings(k):
    if k == 1:
        return ["0","1"]
    res = []
    for x in gen_universal_strings(k-1):
        res.append("0" + x)
        res.append("1" + x)
    return res


def non_branching_path(g):
    res = []
    ig = inverse_graph(g)
    vis = defaultdict(int)
    for v in g.keys():
        if len(g[v]) != 1 or len(ig[v]) != 1:
            vis[v] = 1
            for u in g[v]:
                p = [v,u]
                w = u
                vis[u] = 1
                while len(g[w]) == 1 and len(ig[w]) == 1:
                    w = g[w][0]
                    vis[w] = 1
                    p.append(w)
                res.append(p)
    for v in g.keys():
        if vis[v] != 1:
            w = v
            p = []
            found = False
            while len(g[w]) == 1 and len(ig[w]) == 1:
                vis[w] = 1
                w = g[w][0]
                p.append(w)
                if w == v:
                    found = True
                    break
            if found:
                res.append([v]+p)
    return res