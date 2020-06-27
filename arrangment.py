def greedy_sort(p):
    res = []
    for i in range(1, max(p)+1):
        k=i-1
        if i in p:
            j = p.index(i)
            if j != k:
                p = p[:k] + [-x for x in p[k:j+1]][::-1] + p[j+1:]
                res.append(list(p))
                p[k] = -p[k]
                res.append(list(p))
        else:
            j = p.index(-i)
            if j != k:
                p = p[:k] + [-x for x in p[k:j+1]][::-1] + p[j+1:]
                res.append(list(p))
            else:
                p[k] = -p[k]
                res.append(list(p))

    return res

def breakpoin_count(x):
    x = [0] + x + [len(x)+1]
    return sum([x[i] - x[i-1] != 1 for i in range(1,len(x))])

def chromosome_to_cycle(c):
    res = [0]*2*len(c)
    for j in range(len(c)):
        i = c[j]
        if i > 0:
            res[j*2] = 2*i-1
            res[j*2+1] = 2*i
        else:
            res[j*2+1] = -2*i-1
            res[j*2] = -2*i
    return res

def cycle_to_chromosome(c):
    res=[int(c[j*2+1]/2) if c[j*2] < c[j*2+1] else int(-c[j*2]/2) for j in range(int(len(c)/2))]
    return res            

def colored_edges(p):
    e = []
    for c in p:
        n = chromosome_to_cycle(c)
        for j in range(len(c)):
            e.append((n[j*2-1],n[j*2]))
    return e

def black_edges(p):
    e = []
    if isinstance(p,int):
        n = p
    else:
        n = sum(len(c) for c in p)
    for i in range(1,n+1,2):
        e.append((i, i+1))
        e.append((i+1, i))
    return e

def nodes(g):
    return set(list(g.keys()) + [y for x in g.values()  for y in x])

def edges_to_graph(e, p, g=None, directed=False):
    if not g:
        g = dict()
    for u,v in e:
        if u in g:
            g[u].update({v:p})
        else:
            g[u] = {v:p}
        if not directed:
            if v in g:
                g[v].update({u:p})
            else:
                g[v] = {u:p}
        elif v not in g:
            g[v] = {}
    return g


def merge_graphs(g1, g2):
    g = {}
    for u in g1.keys():
        if u not in g:
            g[u] = {}
        for v in g1[u]:
            g[u][v] = g1[u][v]
    for u in g2.keys():
        if u not in g:
            g[u] = {}
        for v in g2[u]:
            g[u][v] = g2[u][v]
    return g

def graph_to_genome(g):
    vis = set()
    p = []
    for u in g.keys():
        if u in vis:
            continue
        c = []
        has_next = True
        while has_next:
            vis.add(u)
            c.append(u)
            has_next = False
            for v in g[u]:
                if v not in vis:
                    has_next = True
                    u = v
                    break
        if abs(c[0]-c[-1]) == 1:
            c = c[-1:] + c[:-1]
        p.append(cycle_to_chromosome(c))
    return p


def break_distance(p, q):
    g = edges_to_graph(colored_edges(p),1)
    g = edges_to_graph(colored_edges(q),2,g)
    p = graph_to_genome(g)
    n = sum(len(x) for x in p)
    return n - len(p)

def break_graph(gl, i1, i2, i3, i4):
    for i in range(len(gl)):
        x = gl[i]
        if  (x[0], x[1]) == (i1, i2) or (x[0], x[1])  == (i2, i1):
            gl[i] = (i2, i4)
        elif (x[0], x[1])  == (i3, i4) or (x[0], x[1])  == (i4, i3):
            gl[i] = (i1, i3)
    gl = sorted(gl, key= lambda x: x[0])
    return gl

def break_genome(p, i1, i2, i3, i4):
    ce = colored_edges(p)
    be = black_edges(ce)
    g = edges_to_graph(break_graph(ce,i1,i2,i3,i4), 2 , {}, False)
    g = edges_to_graph(be,1, g)
    p = graph_to_genome(g)
    return p

import pprint


def shortest_rearrangement_scenario(p, q):
    r = edges_to_graph(colored_edges(p),2)
    b = edges_to_graph(colored_edges(q),3)
    g = merge_graphs(r, b)
    res = [p]
    for u in g.keys():
        while len(g[u]) != 1:
            i1 = u
            i2 = list(r[u].keys())[0]
            i4 = list(b[u].keys())[0]
            i3 = list(r[i4].keys())[0]
            r[i1] = {i4:2}
            r[i4] = {i1:2}
            r[i2] = {i3:2}
            r[i3] = {i2:2}
            g = merge_graphs(r, b)
            p=break_genome(p,i1,i2,i4,i3)
            res.append(p)
    return res