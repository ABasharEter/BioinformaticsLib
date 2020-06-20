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
    res=[]
    for j in range(int(len(c)/2)):

        if c[j*2] < c[j*2+1]:
            res.append(int(c[j*2+1]/2))
        else:
            res.append(int(-c[j*2]/2))
    return res            

def colored_edges(p):
    e = []
    for c in p:
        n = chromosome_to_cycle(c)
        n = n[1:] + n[:1]
        for j in range(len(c)):
            e.append((n[j*2],n[j*2+1]))
    return e

def nodes(g):
    return set(list(g.keys()) + [y for x in g.values()  for y in x])

def dfs(v, g, r):
    s = [v]
    res = []
    i = 0
    while s:
        v = s.pop()
        if r[v] < 0:
            r[v] = i
            i += 1
            s.append(v)
            for u in g[v]:
                if r[u] < 0:
                    s.append(u)
        else:
            res.append(v)
    return res

def compoents(g):
    r = {k:-1 for k in nodes(g)}
    res = []
    for v in r.keys():
        if r[v] < 0:
            res.append(dfs(v, g, r))
    return res

def edges_to_graph(e, p, g=dict(), directed=False):
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

def graph_to_genom(el):
    g = edges_to_graph(el, {}, True)
    last_node = max(g.keys())
    for i in range(2,last_node+1):
        if len(g[i]) == 0:
            g[i].append(i-1)
    if len(g[1]) == 0:
        g[1].append(last_node)
    print(g)
    cs  = compoents(g)
    p = []
    for c in cs:
        c = set(c)
        print(c)
        c = [y for x in el for y in x if x[0] in c and x[1] ]
        p.append(cycle_to_chromosome(c))
    return p

def break_distance(p, q):
    gp = colored_edges(p)
    gq = colored_edges(q)
    n = max(abs(x) for y in p+q for x in y )
    g = edges_to_graph(gq, edges_to_graph(gp))
    c = compoents(g)
    return n-len(c)

def break_graph(gl, i1, i2, i3, i4):
    for i in range(len(gl)):
        x = gl[i]
        if  (x[0], x[1]) == (i1, i2) or (x[0], x[1])  == (i2, i1):
            gl[i] = (i2, i4)
        elif (x[0], x[1])  == (i3, i4) or (x[0], x[1])  == (i4, i3):
            gl[i] = (i1, i3)
    gl = sorted(gl, key= lambda x: x[0])
    return gl

def break_genom(p, i1, i2, i3, i4):
    g = colored_edges(p)
    #print(g)
    g = break_graph(g, i1, i2, i3, i4)
    #print(g)
    p = graph_to_genom(g)
    return p
    
