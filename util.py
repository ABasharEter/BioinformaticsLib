def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def print_array(a):
    [print("\t".join(list(map(str, map(int, x))))) for x in a]    


def print_weighted_comma_graph(g):
    return "\n".join(["->".join((str(k),",".join(":".join((str(u),str(w))) for u,w in v.items()))) for k,v in g.items()])


def print_weighted_lines_graph(g):
    return "\n".join(["->".join((str(k),":".join((str(u),str(w)))))  for k,v in g.items() for u,w in v.items()])


def graph2edges(g):
    e = []
    for k,v in g.items():
        if isinstance(v, dict):
            v = list(v.items())
        e.extend(zip([k]*len(v),v))
    return e