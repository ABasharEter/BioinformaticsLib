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
    return min(d[i,j] +d[j,k] -d[i,k] for i in range(d.shape[0]) for k in range(i, d.shape[0]) if i!=j and k!=j)/2
