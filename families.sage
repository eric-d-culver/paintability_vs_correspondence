

def family_1(n):
    G = Graph(n)
    cycle = n-3
    edges = []
    for i in [0..cycle-2]:
        edges.append((i,i+1))
    edges.append((cycle-1,0))
    edges.append((0,n-3))
    edges.append((n-3,n-2))
    edges.append((n-2,0))
    for i in [1..n-2]:
        edges.append((n-1,i))
    G.add_edges(edges)
    return G.graph6_string()

def family_2(m):
    G = Graph(2*m + 2)
    edges = []
    for i in [0..m-1]:
        for j in [i+1..m-1]:
            edges.append((i,j))
    for i in [m..2*m-1]:
        for j in [i+1..2*m-1]:
            edges.append((i,j))
    for i in [0..2*m-1]:
        edges.append((i,2*m))
        edges.append((i,2*m+1))
    G.add_edges(edges)
    return G.graph6_string()

#print(family_1(6))
#print(family_1(7))
#print(family_1(8))
#print(family_1(9))
