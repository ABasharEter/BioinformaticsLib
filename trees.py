import numpy as np

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
       