import numpy as np

def coin_excange(money, coins):
    dp = [0]*(money +10)
    for i in range(1, money+3):
        x = [dp[i-c] +1 for c in coins if i-c>=0]
        dp[i] = min(x) if len(x)>0 else np.inf
    return dp[money]

        
def manhaten_tourist(down, right, n, m):
    dp = np.zeros((n+10,m+10))
    for i in range(1,n+1):
        dp[i,0] = dp[i-1,0] + down[i-1,0]
    for i in range(1,m+1):
        dp[0,i] = dp[0,i-1] + right[0,i-1]
    for i in range(1,n+1):
        for j in range(1,m+1):
            dp[i,j] = max(dp[i-1,j] + down[i-1,j], dp[i,j-1] + right[i,j-1])
    return dp[n,m]

    
def pares_array(text):
    return np.array([[float(x) for x in y.split(" ")] for y in text.split("\n")])

def longest_common_subsequnce(s1,s2):
    dp = np.zeros((len(s1)+10,len(s2)+10))
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            dp[i,j] = max(dp[i-1,j], dp[i,j-1], dp[i-1,j-1] + (s1[i-1] == s2[j-1]))
    res = []
    i = len(s1)
    j = len(s2)
    while i>0 and j>0:
        if s1[i-1] == s2[j-1]:
            res.append(s1[i-1])
            i-=1
            j-=1
        elif dp[i-1,j] < dp[i,j-1]:
            j-=1
        else:
            i-=1
    return "".join(reversed(res))

def longest_common_subsequnce3(s1,s2,s3):
    dp = np.zeros((len(s1)+1,len(s2)+1,len(s3)+1))
    for i in range(len(s1)+1):
        for j in range(len(s2)+1):
            for k in range(len(s3)+1):
                if i==0 and j==0 and k==0:
                    continue
                elif (s1[i-1] == s2[j-1] and s1[i-1] == s3[k-1]):
                    dp[i,j,k] = dp[i-1,j-1,k-1]+1
                else:
                    dp[i,j,k] = max(dp[i-1,j,k], dp[i,j-1,k], dp[i,j,k-1])
    res = []
    i = len(s1)
    j = len(s2)
    k = len(s3)
    print(dp[i,j,k])
    while i>0 and j>0 and k>0:
        if dp[i,j,k] == dp[i-1,j-1,k-1]+1:
            res.append(s1[i-1])
            i-=1
            j-=1
            k-=1
        if dp[i,j,k] == dp[i-1,j,k]:
            i -=1
        elif dp[i,j,k] == dp[i,j-1,k]:
            j -=1
        elif  dp[i,j,k] == dp[i,j,k-1]:
            k -=1
    return "".join(reversed(res))



def topological_sort(dag):
    nodes = list(set(list(dag.keys()) + [x for v in dag.values() for x in v.keys()]))
    q = []
    in_deg = {u:0 for u in nodes}
    for u,adj in dag.items():
        for v in adj:
            in_deg[v] +=1
    top_sort = []
    for u in nodes:
        if in_deg[u] == 0:
            q.append(u)
    while q:
        u = q.pop()
        top_sort.append(u)
        for v in dag[u]:
            in_deg[v] -=1
            if in_deg[v] == 0:
                q.append(v)
    return top_sort

def parse_weighted_edige_list(text):
    edges = [x for x in text.split("\n") if len(x) > 0]
    g = {}
    for e in edges:
        u,l = [x.strip() for x in e.split("->")]
        v,w = [x.strip() for x in l.split(":")]
        if v not in g:
            g[v] = {}
        if u not in g:
            g[u] = {}
        
        g[u][v] = int(w)
    return g

def longest_path(s,e,dag):
    nodes = list(set(list(dag.keys()) + [x for v in dag.values() for x in v.keys()]))
    top_sort = topological_sort(dag)
    if len(nodes) != len(top_sort):
        return np.inf
    dp = {k:0 for k in top_sort}
    bp = {k:None for k in top_sort}
    for v in top_sort:
        for u,w in dag[v].items():
            if dp[u] < dp[v]+w:
                bp[u] = v
                dp[u] = dp[v]+w
    if not e:
        e = [x for x in dp.keys() if dp[x] == max(dp.values())][0]
    p = [e]
    while bp[p[-1]] and p[-1] != s:
        p.append(bp[p[-1]])
    if not s:
        s = p[-1]
    return dp[e] - dp[s], list(reversed(p))
    
def parse_scoreing_matrix(m):
    m = [[v for v in l.split(" ") if len(v)> 0] for l in m.split("\n")]
    c = m[0]
    w = {u:{v:0 for v in c} for u in c}
    for i in range(1,len(m)):
        for j in range(1,len(m[i])):
            w[m[i][0]][m[0][j-1]] = int(m[i][j])
    return w

def constant_edit_score(i_e,d_e,i_o,d_o,c,m):
    res = {a:{b:c if a != b else m for b in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"} for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
    res['-'] = {a:i_o for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
    res['*'] = {a:i_e for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
    for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        res[a]['-'] = d_o
        res[a]['*'] = d_e
    return res

def add_score(w, char, insersion_score, deletion_score):
    w[char] = {k:insersion_score for k in w.keys()}
    for k in w.keys():
        w[k][char] = deletion_score

def matching_alignment(s1, s2, w, alignment_type='global'):
    # - is sigma, * is epsilon
    dp = np.zeros((len(s1)+1,len(s2)+1))
    dp_lower = np.zeros((len(s1)+1,len(s2)+1))
    dp_upper = np.zeros((len(s1)+1,len(s2)+1))
    for i in range(1,len(s1)+1):
        dp[i,0] = w[s1[i-1]]["-"] +  w[s1[i-1]]["*"]*(i-1)
        if alignment_type=='local' or alignment_type=='fitting' or alignment_type=='overlap':
            dp[i,0] = max(dp[i,0],0)
        dp_lower[i,0] = dp[i][0]
        dp_upper[i,0] = -np.inf

    for j in range(1,len(s2)+1):
        dp[0,j] = w["-"][s2[j-1]] + w["*"][s2[j-1]]*(j-1)
        if alignment_type=='local':
            dp[0,j] = max(dp[0,j],0) 
        dp_lower[0,j] = -np.inf
        dp_upper[0,j] = dp[0,j]
        
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            dp_lower[i,j] = max(dp_lower[i-1,j] + w[s1[i-1]]["*"], dp[i-1,j] + w[s1[i-1]]["-"])
            dp_upper[i,j] = max(dp_upper[i,j-1] + w["*"][s2[j-1]],  dp[i,j-1] + w["-"][s2[j-1]])
            dp[i,j] = max(dp_upper[i,j], dp_lower[i,j], dp[i-1,j-1] + (w[s1[i-1]][s2[j-1]]))

            if alignment_type=='local':
                dp[i,j] = max(0,dp[i,j])
    layer = 'middle'
    if alignment_type=='local':
        i,j = [(i,j) for i in range(len(s1)+1) for j in range(len(s2)+1) if dp[i,j] == m][0]
        m = np.max(dp)
    elif alignment_type == 'fitting':
        j = len(s2)
        m = np.max(dp[:len(s1)+1,j])
        i = [i for i in range(len(s1)+1) if dp[i,j] == m][0]
    elif alignment_type == 'overlap':
        i = len(s1)
        m = np.max(dp[i,:len(s2)+1])
        j = [j for j in range(len(s2)+1) if dp[i,j] == m][0]
    elif alignment_type=='global':
        i = len(s1)
        j = len(s2)
        m = dp[len(s1),len(s2)]
    r1,r2 = [],[]
    while (i>0 or j>0):
        if layer == 'middle':
            if i>0 and j>0 and dp[i,j] == dp[i-1,j-1] + (w[s1[i-1]][s2[j-1]]):
                r1.append(s1[i-1])
                r2.append(s2[j-1])
                i -= 1
                j -= 1
            elif dp[i,j] == dp_lower[i,j]:
                layer = 'lower'
            else:
                layer = 'upper'
        elif layer == 'lower':
            if i>0:
                r1.append(s1[i-1])
                r2.append('-')
                i -= 1
                if dp_lower[i,j] == dp[i-1,j] + w[s1[i-1]]["-"]:
                    layer = 'middle'
            else:
                layer = 'middle'
        elif layer == 'upper':
            if j>0:
                r1.append('-')
                r2.append(s2[j-1])
                j -= 1
                if dp_upper[i,j] == dp[i,j-1] + w["-"][s2[j-1]]:
                    layer = 'middle'
            else:
                layer = 'middle'
        if dp[i,j] == 0 and alignment_type == 'local':
            break
        if j == 0 and (alignment_type == 'fitting' or alignment_type == 'overlap'):
            break
        
    return m ,"".join(reversed(r1)),"".join(reversed(r2))


def middle_edge(s1, s2, w):
    # - is sigma, * is epsilon
    if len(s1) <= 1 or len(s2) <= 1:
        s,r1,r2 = matching_alignment(s1, s2, w, 'global')
    else:
        n = len(s1)
        n2 = int(n/2)
        score_l = np.array(_matching_alignment_score(s1[:n2], s2, w))
        score_r = np.array(_matching_alignment_score("".join(reversed(s1[n2:])), "".join(reversed(s2)), w)[::-1])
        score = score_l + score_r
        m2 = np.argmax(score)
        s = score[m2]
        _,r1_1,r2_1 = matching_alignment2(s1[:n2], s2[:m2], w)
        _,r1_2,r2_2 = matching_alignment2(s1[n2:], s2[m2:], w)
        r1,r2 = r1_1+r1_2, r2_1+r2_2
    return s,r1,r2

def matching_alignment2(s1, s2, w):
    # - is sigma, * is epsilon
    if len(s1) <= 1 or len(s2) <= 1:
        s,r1,r2 = matching_alignment(s1, s2, w, 'global')
    else:
        n = len(s1)
        n2 = int(n/2)
        score_l = np.array(_matching_alignment_score(s1[:n2], s2, w))
        score_r = np.array(_matching_alignment_score("".join(reversed(s1[n2:])), "".join(reversed(s2)), w)[::-1])
        score = score_l + score_r
        m2 = np.argmax(score)
        s = score[m2]
        _,r1_1,r2_1 = matching_alignment2(s1[:n2], s2[:m2], w)
        _,r1_2,r2_2 = matching_alignment2(s1[n2:], s2[m2:], w)
        r1,r2 = r1_1+r1_2, r2_1+r2_2
    return s,r1,r2

def _matching_alignment_score(s1, s2, w):
    dp = np.zeros((2,len(s2)+1))
    dp_lower = np.zeros((2,len(s2)+1))
    dp_upper = np.zeros((2,len(s2)+1))
    for j in range(1,len(s2)+1):
        dp[0,j] = w["-"][s2[j-1]] + w["*"][s2[j-1]]*(j-1)
        dp_lower[0,j] = -np.inf
        dp_upper[0,j] = dp[0,j]
    layer = 0
    for i in range(1,len(s1)+1):
        layer = 1-layer
        dp[layer,0] = w[s1[i-1]]["-"] +  w[s1[i-1]]["*"]*(i-1)
        dp_lower[layer,0] = dp[layer,0]
        dp_upper[layer,0] = -np.inf
        for j in range(1,len(s2)+1):
            dp_lower[layer,j] = max(dp_lower[1-layer,j] + w[s1[i-1]]["*"], dp[1-layer,j] + w[s1[i-1]]["-"])
            dp_upper[layer,j] = max(dp_upper[layer,j-1] + w["*"][s2[j-1]],  dp[layer,j-1] + w["-"][s2[j-1]])
            dp[layer,j] = max(dp_upper[layer,j], dp_lower[layer,j], dp[1-layer,j-1] + (w[s1[i-1]][s2[j-1]]))
    return dp[layer]