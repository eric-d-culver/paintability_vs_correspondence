
n = 6

graph6_strings = []

# generate the interesting graphs
for G in graphs(n):
    if G.is_connected() and min(G.degree_iterator()) >= 3:
	Delta = max(G.degree_iterator())
        if Delta > 3:
	    omega = G.clique_number()
            if omega < Delta:
                #graph6_strings.append(G.graph6_string())
                i = G.chromatic_number()
                print G.graph6_string(), i

# extract the subgraph minimal graphs of this set 
# NOTE: subgraph_search in Sage is slow, so disable if too slow
#minimals = []
#for graph6 in graph6_strings:
    #G = Graph(graph6)
    #new_minimal = True
    #for min_graph6 in minimals:
        #if G.subgraph_search(Graph(min_graph6)) is not None:
            #new_minimal = False
    #if new_minimal:
        #minimals.append(graph6)
        #print(graph6)
