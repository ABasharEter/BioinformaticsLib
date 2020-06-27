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
    txt = """(+1 +2 +3 +4)(+5 +6)(+7+8+9)
""".split("\n")
    r = pares_input(txt[0])
    #b = pares_input(txt[1])
    z = chromosome_to_cycle(r,b)
    print(z)
    #print(", ".join((f"({u}, {v})" for u,v in y)))
    y = [" ".join(["(" + " ".join([str(x) if x < 0 else "+" + str(x) for x in v]) + ")" for v in yy]) for yy in z]
    print("\n".join(y))
    
    
    