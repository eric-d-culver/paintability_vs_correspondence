
n = 8

graph6_strings = []

# generate the interesting graphs
for G in graphs(n):
    if G.is_connected() and min(G.degree_iterator()) >= 3:
	Delta = max(G.degree_iterator())
        c = max(G.cores()) + 1 # coloring number of graph
        ub = min(c, Delta) # c <= Delta always
        if ub > 3: 
	    omega = G.clique_number()
            mad = G.maximum_average_degree()
            lb = max(omega, ceil(mad/2)) # AT(G) > mad(G)/2
            if ub - lb > 1: # want gap of at least 2
                #graph6_strings.append(G.graph6_string())
                i = G.chromatic_number()
                lb = max(i, lb)
                if ub - lb > 1: # want gap of at least 2
                    print G.graph6_string(), lb, i, ub

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
