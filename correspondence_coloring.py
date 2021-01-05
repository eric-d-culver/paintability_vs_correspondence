#!/usr/bin/python3

# basic idea:
# represent a graph as a spanning tree plus a bunch of nontree edges.
# assume the correspondences on the spanning tree edges are all straight. 
# (reduces the number of correspondences to consider by a sizeable amount)
# Let k be the number of colors.
# We consider all the possible correspondences on the nontree edges (loop through permutations of k)
# For each, we loop through possible colorings of the verts incident to the nontree edges which respect those correspondences.
# Therefore, we have a tree with some precolored vertices, and we want to know if we can extend that coloring to the entire tree.
# We do this through a bottom-up approach. (the buzzword is "dynamic programming".)
# we build up from the leaves to the root of our tree the number of possible valid colorings with that precoloring.
# At each node, we have a vector of k numbers. 
# The ith element of the vector is the number of valid colorings of the subtree rooted at that node where that node is colored i.
# Rules for building up these vectors:
# Leaf not precolored: vector of all ones
# Leaf precolored i: vector of all zeroes except in the ith position where we have a one.
# Branch not precolored: the ith element is the sum over all colorings of the children which don't include i 
# of the product of all the ways to color the children that coloring.
# Branch precolored i: same as above except we zero out all except the ith element.
# When we reach the root, the sum of the elements of the root's vector gives the number of valid colorings
# of the tree with that precoloring.

import itertools as it

#represent tree as a tuple (node_number, tuple_of_subtrees)

# precolored vertices given by a dictionary with keys of node_numbers and values of colors

def colorings_vector(k, tree, precolored_dict):
    if tree[1] is ():
        # leaf
        if tree[0] in precolored_dict:
            return tuple(int(i == precolored_dict[tree[0]]) for i in range(k))
        else:
            return tuple(1 for i in range(k))
    else:
        # branch
        num_children = len(tree[1])
        # collect the results for the child subtrees
        children = []
        for subtree in tree[1]:
            children.append(colorings_vector(k, subtree, precolored_dict))
        vals = []
        for i in range(k):
            ways = 0
            for colors in it.product(range(k), repeat=num_children):
                if i in colors:
                    continue
                else:
                    prod = 1
                    for j in range(num_children):
                        prod *= children[j][colors[j]]
                    ways += prod
            vals.append(ways)
        if tree[0] in precolored_dict:
            # in this case we wasted time computing the other elements of vals, but I think it is not that big a deal
            return tuple((vals[i] if i == precolored_dict[tree[1]] else 0) for i in range(k))
        else:
            return tuple(vals)

def precolored_tree(k, tree, precolored):
    return sum(colorings_vector(k, tree, precolored))

# represent nontree edges as tuples (node_number, node_number)

# represent correspondences as k-tuples one for each edge.
# edge (u,v) has a correspondence (2,3,1) means that 
# color 1 on u is matched with color 2 on v,
# color 2 on u is matched with color 3 on v, and
# color 3 on u is matched with color 1 on v.
# direction of tuple for edge gives direction of permutation for correspondence.

def valid_precolorings(k, nontree_edges, correspondences):
    precolored_verts = list(set(e[0] for e in nontree_edges).union(set(e[1] for e in nontree_edges)))
    num_precolored = len(precolored_verts)
    for coloring in it.product(range(k), repeat=num_precolored):
        coloring_dict = {key : value for (key, value) in zip(precolored_verts, coloring)}
        for i in range(num_precolored):
            for j in range(len(nontree_edges)):
                if nontree_edges[j][0] != precolored_verts[i]:
                    continue
                else:
                    if correspondences[j][colorings[i]] == colorings[nontree_edges[j][1]]:
                        pass #bad coloring
            pass
    pass

def num_colorings(k, tree, nontree_edges):
    colorings = 0
    for correspondences in it.product(it.permutations(k), repeat=len(nontree_edges)):
        for precoloring in valid_precolorings(k, nontree_edges, correspondences):
            colorings += precolored_tree(k, tree, precoloring)
    return colorings

print(colorings_vector(2, (1, ((2, ()), (3, ((4, ()),)))), {2: 0, 4: 1}))
print(colorings_vector(2, (1, ((2, ()), (3, ((4, ()),)))), {2: 1, 4: 1}))
