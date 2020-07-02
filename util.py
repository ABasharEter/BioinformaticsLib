def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def print_array(a):
    [print("\t".join(list(map(str, map(int, x))))) for x in a]
    