from arrangment import *
def pares_input(text):
    p = text.strip("()").split(")(")
    p = [[int(x) for x in y.split(" ")] for y in p]
    return p

def pares_input2(text):
    p = text.split("), ")
    p = [[int(x.strip("() ")) for x in y.split(",")] for y in p]
    return p

if __name__ == "__main__":
    txt = """(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)
""".split("\n")
    p = pares_input2(txt[0])
    #i1,i2,i3,i4 = [int(x) for x in txt[1].split(", ")]
    y = graph_to_genom(p)
    #print(", ".join((f"({u}, {v})" for u,v in y)))
    y = ["("+" ".join(str(x) if x < 0 else f"+{x}" for x in v)+")" for v in y]
    print(" ".join(y))
    
    
    