import kmer

import numpy as np
import copy
import amino

def parse_wieghtd_graph(data):
    res = {}
    for line in data.split("\n"):
        k,v = line.split("->")
        v = v.split(",")
        k = k.strip()
        res[k] = [x.strip() for x in v]
        for x in v:
            u,w = x.strip().split
            if u not in res:
                res[u] = {}
    return res


def dfs(p, u, g, d, ids):
    for v,x in g[u].items():
        if p and v == p:
            continue
        d[ids[v]] = d[ids[u]] + x
        dfs(u, v, g, d, ids)

def tree2distance(g):
    d = np.zeros((len(g),len(g)))
    keys_ids = {u:i for i,u in enumerate(sorted(g.keys())) }
    leafs_ids = [i for u,i in keys_ids.items() if len(g[u]) == 1]
    for u in g.keys():
        dfs(None, u, g, d[keys_ids[u]], keys_ids)
    return d[leafs_ids,:][:,leafs_ids]


def limb_length(j,d):
    v = [(int((d[i,j] +d[j,k] -d[i,k])/2),i,k) for i in range(d.shape[0]) for k in range(i, d.shape[0]) if i!=j and k!=j]
    m = min(x[0] for x in v)
    r = [x for x in v if x[0] == m][0]
    return r

def addNode(t, j, limb_len, i, k, x):
    l = len(t)
    dist = [float('inf')] * l
    parent = [-1] * l
    q = []
    dist[i] = 0
    q.append(i)
    while len(q) > 0:
        u = q.pop(0)
        for v, w in t[u].items():
            if float('inf') == dist[v]:
                dist[v] = dist[u] + w
                parent[v] = u
                q.append(v)
                if v == k:
                    p = v
                    while x < dist[p]:
                        c = p
                        p = parent[p]
                    if x == dist[p]:
                        t[p][j] = limb_len
                        t[j][p] = limb_len
                    else:
                        nu = max(t.keys()) + 1
                        t[nu] = {j:limb_len}
                        t[j][nu] = limb_len
                        del t[p][c]
                        del t[c][p]
                        t[p][nu] = x-dist[p]
                        t[nu][p] = x-dist[p]
                        t[c][nu] = dist[c]-x
                        t[nu][c] = dist[c]-x
                    return

def additive_phylogeny(d):
    n = d.shape[0]
    t = {k:dict() for k in range(n)}
    t[0][1] = d[0][1]
    t[1][0] = d[1][0]
    for j in range(2, n):
        l, i, k = limb_length(j, d[:j+1,:j+1])
        x = d[i][j] - l
        addNode(t, j, l, i, k, x)
    return t


def delete_ij(D,i,j,cid,cz,k):
    n = D.shape[0]
    if i>j:
        i,j=j,i
    cz[k]= cz[cid[i]]+cz[cid[j]]
    new_D = D
    new_D = np.delete(new_D, (i,j), axis=0)
    new_D = np.delete(new_D, (i,j), axis=1)
    new_row = (D[i,:]*cz[cid[i]]+D[j,:]*cz[cid[j]])/(cz[k])
    new_row[:n-2] = np.delete(new_row, (i,j))
    new_row[n-2] = 0
    new_D = np.vstack((new_D, new_row[:n-2]))
    new_D = np.hstack((new_D, new_row[:n-1].reshape(-1,1)))
    for x in range(n-1):
        if x < j:
            cid[x] = cid[x]
        else:
            cid[x] = cid[x+1]
    for x in range(n-1):
        if x < i:
            cid[x] = cid[x]
        else:
            cid[x] = cid[x+1]
    cid[n-2]=k
    return new_D

def UPGMA(D):
    n = D.shape[0]
    g = {i:{} for i in range(n)}
    age = {i:0 for i in range(n)}
    cid = {i:i for i in range(n)}
    cz = {i:1 for i in range(n)}
    for _ in range(n-1):
        i,j = np.unravel_index((D + np.identity(n)*1e9).argmin(), D.shape)
        d = D[i,j]
        k = len(g)
        age[k] = d/2
        g.setdefault(cid[i],{}).update({k:age[k]-age[cid[i]]})
        g.setdefault(k,{}).update({cid[i]:age[k]-age[cid[i]]})
        g.setdefault(cid[j],{}).update({k:age[k]-age[cid[j]]})
        g.setdefault(k,{}).update({cid[j]:age[k]-age[cid[j]]})
        D = delete_ij(D,i,j,cid,cz,k)
        n = D.shape[0]
    return g,age


def D_star(D):
    t = np.repeat(np.sum(D,axis=0, keepdims=True),D.shape[1], axis=0)
    D = (D.shape[0] -2)*D - t - t.T
    np.fill_diagonal(D,0)
    return D
    

def delete_ij2(D,i,j,cid,k):
    n = D.shape[0]
    if i>j:
        i,j=j,i
    new_D = D
    new_D = np.delete(new_D, (i,j), axis=0)
    new_D = np.delete(new_D, (i,j), axis=1)
    new_row = np.array([D[i,k]+D[j,k] - D[i,j] for k in range(D.shape[0])])
    new_row = new_row/2
    new_row[:n-2] = np.delete(new_row, (i,j))
    new_row[n-2] = 0
    new_D = np.vstack((new_D, new_row[:n-2]))
    new_D = np.hstack((new_D, new_row[:n-1].reshape(-1,1)))
    for x in range(n-1):
        if x < j:
            cid[x] = cid[x]
        else:
            cid[x] = cid[x+1]
    for x in range(n-1):
        if x < i:
            cid[x] = cid[x]
        else:
            cid[x] = cid[x+1]
    cid[n-2]=k
    return new_D

def neighbor_joining(D):
    n = D.shape[0]
    g = {i:{} for i in range(n)}
    cid = {i:i for i in range(n)}
    for _ in range(n-1):
        Ds = D_star(D)
        i,j = np.unravel_index((Ds + np.identity(n)*1e9).argmin(), Ds.shape)
        d = Ds[i,j]
        k = len(g)
        if n > 2:
            td = np.sum(D, axis=0)
            delta = (td[i] - td[j])/(n-2)
            lli = (D[i,j] + delta)/2
            llj = (D[i,j] - delta)/2
            g.setdefault(cid[i],{}).update({k:lli})
            g.setdefault(k,{}).update({cid[i]:lli})
            g.setdefault(cid[j],{}).update({k:llj})
            g.setdefault(k,{}).update({cid[j]:llj})
            D = delete_ij2(D,i,j,cid,k)
            n = D.shape[0]
        else:
            lli = llj = D[0,1]
            g.setdefault(cid[i],{}).update({cid[j]:D[0,1]})
            g.setdefault(cid[j],{}).update({cid[i]:D[0,1]})
    return g

def read_bi_parsimony_tree_array(txt, rooted =False):
    txt = txt.split("\n")
    g = {}
    s = set()
    for l in txt[1:]:
        u, v = l.strip().split("->")
        g.setdefault(u,set()).add(v)
        g.setdefault(v,set()).add(u)
        if u.isalpha():
            s.add(u)
        if v.isalpha():
            s.add(v)
    n = int(txt[0])
    q = list(s)
    r = None
    vis = set()
    while q:
        u = q.pop(0)
        if u in vis:
            continue
        vis.add(u)
        r = u
        q.extend(g[u])
    if not rooted:
        u = None
        for v in g[r]:
            if not v.isalpha():
                u = v
                break
        g[r].remove(u)
        g[u].remove(r)
        nr = '#'
        g[nr]= [r,u]
        r = nr
    q = [r]
    t = {}
    vis = set()
    idxs = {r:1}
    while q:
        u = q.pop(0)
        if u in vis:
            continue
        vis.add(u)
        if u.isalpha():
            t[idxs[u]] = u
        else:
            print(u,g[u])
            x,y = [i for i in g[u] if i not in vis]
            idxs[x] = idxs[u]*2
            idxs[y] = idxs[u]*2+1
            t[idxs[u]] = None
        q.extend(g[u])
    return t
very_larg = 1e18
from pprint import pprint
def small_parsimony(ta):
    n = len([x for x in ta.values() if x][0])
    tao = {k:False if v else True for k,v in ta.items()}
    ta = {k:list(v) if v else [None]*n for k,v in ta.items()}
    for i in range(n):
        tp = {k:{} for k in ta.keys()}
        for j in reversed(sorted(ta.keys())):
            if tao[j]:
                tp[j] = {k:min([v + (0 if k2 == k else 1) for k2,v in tp[j*2+1].items()]) + min([v + (0 if k2 == k else 1) for k2,v in tp[j*2].items()])
                                 for k in kmer.dna_nucleotides}
            else:
                tp[j] = {k:(very_larg if k != ta[j][i] else 0) for k in kmer.dna_nucleotides}
        #if i in (6,7,8):
            #print(i)
            #pprint(tp)
        for j in sorted(ta.keys()):
            m = min(tp[j].items(), key=lambda x: x[1])[0]
            if int(j/2) > 0 and tp[j][m] >= tp[j][ta[int(j/2)][i]]:
                ta[j][i] = ta[int(j/2)][i]
            else:
                ta[j][i] = m
    ta = {k:''.join(v) for k,v in ta.items()}
    return ta

def print_ta(ta,rooted=False):
    s = 0
    txt = ""
    for i in ta.keys():
        if i*2 not in ta or i*2+1 not in ta:
            continue
        if not rooted and i == 1:
            continue
        d = kmer.hamming(ta[i],ta[i*2])
        s += d
        txt += f"{ta[i]}->{ta[i*2]}:{d}\n"
        txt += f"{ta[i*2]}->{ta[i]}:{d}\n"
        d = kmer.hamming(ta[i],ta[i*2+1])
        s += d
        txt += f"{ta[i]}->{ta[i*2+1]}:{d}\n"
        txt += f"{ta[i*2+1]}->{ta[i]}:{d}\n"
    if not rooted:
        d = kmer.hamming(ta[2],ta[3])
        s += d
        txt += f"{ta[2]}->{ta[3]}:{d}\n"
        txt += f"{ta[3]}->{ta[2]}:{d}\n"
    return s,txt

def ta2g(ta,rooted=False):
    ta = [None] + ta
    g = {}
    s = 0
    for i in range(1,int(len(ta)/2)):
        g.setdefault(ta[i],{}).update({
            ta[i*2]:kmer.hamming(ta[i],ta[i*2]),
            ta[i*2+1]:kmer.hamming(ta[i],ta[i*2+1]),
        })
        if i >0 or rooted:
            s += g[ta[i]][ta[i*2]] + g[ta[i]][ta[i*2+1]]
        g.setdefault(ta[i*2], {}).update({ta[i]:g[ta[i]][ta[i*2]]})
        g.setdefault(ta[i*2+1], {}).update({ta[i]:g[ta[i]][ta[i*2+1]]})
    if not rooted:
        del g[ta[1]]
        del g[ta[2]][ta[1]]
        del g[ta[3]][ta[1]]
        g[ta[2]][ta[3]] = g[ta[3]][ta[2]] = kmer.hamming(v,u)
        s += g[ta[2]][ta[3]]
    return s,g

def nearest_neighbors(t,e):
    u,v = e
    t1 = copy.deepcopy(t)
    t2 = copy.deepcopy(t)
    xu,yu = [i for i in t[u] if i != v]
    xv,yv = [i for i in t[v] if i != u]
    t1[u] = [xu,xv,v]
    t1[v] = [yu,yv,u]
    t1[yu].remove(u)
    t1[yu].append(v)
    t1[xv].remove(v)
    t1[xv].append(u)
    t2[u] = [xu,yv,v]
    t2[v] = [xv,yu,u]
    t2[yv].remove(v)
    t2[yv].append(u)
    t2[yu].remove(u)
    t2[yu].append(v)
    return t1,t2
    
    
    
def paths_dfs(g, path, paths = []):   
    u,w = path[-1]              
    if u in g and g[u]:
        for k,v in g[u].items():
            new_path = path + [(k,v)]
            paths = paths_dfs(g, new_path, paths)
    else:
        paths += [path]
    return paths

def spec2g(s):
    imassmap = {v:k for k,v in amino.amino_mass.items()}
    s = [0] + s
    g = {k:{} for k in s}
    for i in range(1,len(g)):
        for j in range(i):
            d = s[i] - s[j]
            if d in imassmap:
                g[s[j]][s[i]] = imassmap[d]
    return g

def ideal_spec(p):
    s = [amino.amino_mass[x] for x in p]
    r = []
    for i in range(1,len(s)):
        r.append(sum(s[:i]))
        r.append(sum(s[i:]))
    r.append(sum(s))
    return sorted(r)
        
    

def decoding_ideal_spectrum(s):
    g = spec2g(s)
    p = paths_dfs(g,[(0,None)])
    for pi in p:
        pi = [x[1] for x in pi[1:]]
        if ideal_spec(pi) == sorted(s):
            return "".join(pi)

def pep2vec(p):
    s = [amino.amino_mass[x] for x in p]
    s = [sum(s[:i]) for i in range(1,len(s)+1)]
    a = [0]*(max(s)+1)
    for i in s:
        a[i] = 1
    return a[1:]

def vec2pep(v):
    s = [0] + [i+1 for i in range(len(v)) if v[i] ==1]
    imassmap = {v:k for k,v in amino.amino_mass.items()}
    p = "".join([imassmap[s[i]-s[i-1]] for i in range(1,len(s))])
    return p
    
