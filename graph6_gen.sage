
n = 6

for G in graphs(n):
    if G.is_connected() and min(G.degree_iterator()) >= 3:
	Delta = max(G.degree_iterator())
	omega = G.clique_number()
        if omega < Delta:
            print(G.graph6_string())
