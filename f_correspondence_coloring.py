#!/usr/bin/python3

# Represent graph as list of edges, which are pairs of numbers (want to include number of verts?)

# Represent partial correspondence as list of lists (or tuples), where the ith element is what i is mapped to. Since we fill in the perms from 1 to n, this will work.

import sys
from functools import lru_cache
from itertools import count, filterfalse
import math

from nauty_utils import graph6_file, graph6_to_dict


# tree is a tuple, first element is label of root, second element is list of subtrees
# empty tree is empty tuple

# helper function
@lru_cache(maxsize=None)
def tree_neighbors(v, tree):
    res = []
    if tree[0] == v: # v is root
        for u, _ in tree[1]:
            res.append(u)
        return res
    for u, subsubtrees in tree[1]:
        if u == v: # v is root of immediate subtree
            res = [tree[0]]
            for w, _ in subsubtrees:
                res.append(w)
            return res
    # else v is possibly lower in the tree
    for x in tree[1]:
        res = tree_neighbors(v, x)
        if len(res) != 0:
            return res
    return [] # no subtrees

# helper function
@lru_cache(maxsize=None)
def neighbors(v, tree, nontree):
    res = []
    for x in tree_neighbors(v, tree):
        res.append((x, 0))
    for x,y in nontree:
        if x == v:
            res.append((y, -1))
        elif y == v:
            res.append((x, 1))
    return res

# helper function
@lru_cache(maxsize=None)
def tree_vertices(tree):
    return [tree[0]] + sum(map(tree_vertices, tree[1]), [])

# helper function
@lru_cache(maxsize=None)
def vertices(tree, nontree):
    res = []
    res = res + tree_vertices(tree)
    for x,y in nontree:
        if x not in res:
            res.append(x)
        if y not in res:
            res.append(y)
    return res

# helper function
@lru_cache(maxsize=None)
def next_missing(A, a):
    return next(filterfalse(set(A).__contains__, count(a+1)))

# takes partial correspondence and yields partial correspondence with edge added to it
# should loop through and yield all the possible ways to add an edge
def add_correspondence_loop(perm, colors):
    for i in range(colors):
        if i not in perm:
            yield perm + (i,)

def everything_but(i, k):
    return list(j for j in range(k) if j != i)

# c_uv is color of given vertex
# orient is whether the perm goes away (+1) or towards (-1) the given vertex
# perm is partial correspondence on edge
# returns list of color choices for neighboring vertex
def color_choices(c_uv, orient, perm, k):
    if orient == 1:
        if c_uv < len(perm):
            return everything_but(perm[c_uv], k) # everything except what c_uv corresponds to
        else: # things that correspond to something
            return list(perm)
    elif orient == -1:
        if c_uv in perm:
            return everything_but(perm.index(c_uv), k) # everything except what c_uv corresponds to
        else: # things that correspond to something
            return list(range(len(perm)))
    elif orient == 0:
        return everything_but(c_uv, k) # everything but c_uv

# c_u, c_v, colors of endpoints
# orient does perm go uv or vu
# perm is partial correspondence on edge
# returns boolean of whether this coloring works for edge
def color_test(c_u, c_v, orient, perm):
    if orient == 1:
        if c_u < len(perm): # c_u corresponds to something
            return perm[c_u] != c_v # so it can't be c_v
        else: # c_u does not correspond to anything
            return c_v in perm # so c_v should correspond to something
    elif orient == -1:
        if c_v < len(perm): # c_v corresponds to something
            return perm[c_v] != c_u # so it can't be c_u
        else: # c_v does not correspond to anything 
            return c_u in perm # so c_u should correspond to something
    elif orient == 0:
        return c_u != c_v


# takes graph, partial correspondence, partial coloring, and list of uncolored vertex and attempts to extend coloring to the uncolored vertices
# tree: tuple of root and list of subtrees (empty tree is empty tuple)
# nontree: list of edges as pairs
# partial correspondence: dict of tuples, keys are edges, values are tuples giving partial correspondence for that edge
# partial coloring: dict of number, keys are vertices (also numbers), values are color of that vertex
# uncolored: list of uncolored vertices
# colors: dictionary giving number of colors available at each vertex
# returns extended coloring if exists, else returns None
# algorithm works by attempting to extend coloring to first element of list, then recursively calling itself on rest of list
def extend_coloring(tree, nontree, correspondences, coloring, uncolored, colors):
    if len(uncolored) == 0: 
        # if everything is colored, we are done
        return coloring
    # otherwise, attempt to extend coloring to first element of list, then recurse
    new_coloring = dict(coloring)
    vertex = uncolored[0]
    remaining = uncolored[1:]
    for i in range(colors[vertex]):
        # test if i is good color for vertex
        good_color = True
        for u, orient in neighbors(vertex, tree, nontree):
            if u in coloring:
                if orient == 1:
                    if not color_test(coloring[u], i, orient, correspondences[(u, vertex)]):
                        good_color = False
                        break
                elif orient == -1:
                    if not color_test(i, coloring[u], orient, correspondences[(vertex, u)]):
                        good_color = False
                        break
                elif orient == 0:
                    if not color_test(coloring[u], i, orient, []):
                        good_color = False
                        break
        if good_color:
            # i is good color for vertex
            new_coloring[vertex] = i
            # recursively color rest of graph
            extended_coloring = extend_coloring(tree, nontree, correspondences, new_coloring, remaining, colors)
            if extended_coloring is not None:
                return extended_coloring
    return None

# bad correspondence algorithm:
# first, set all edges to partial correspondence (0,)
# then loop:
# attempt to color
# if successful, step_correspondence on current edge 
# if unsuccessful, add_correspondence on current edge
# move forward an edge (to distribute moves over all edges)
# how does back tracking work?

# function takes graph, partial correspondence, and list of edges
# try to color graph with partial correspondence
# if successful, return None
# if unsuccessful: 
# if list of edges is empty, return correspondence (it is bad correspondence)
# else, first edge in list is current edge
# add_correspondence on current edge
# if current edge is not full, then move to end of list
# recursively call function on graph with new correspondences and list
# if returns something, return it
# if returns None, step_correspondence on current edge and recurse again (so really we want a loop, meaning we need a list of the possible add_correspondences)
# if they all return None, return None yourself

def bad_correspondence(tree, nontree, correspondences, edges, colors, depth, count):
    #print("Considering correspondence", correspondences)
    # shortcut if some edge of correspondences is ()
    if all(partial != () for partial in correspondences):
        test_coloring = extend_coloring(tree, nontree, correspondences, dict(), vertices(tree, nontree), colors)
        if test_coloring is not None:
            if verbosity > 1:
                print("Coloring", test_coloring, "for correspondences", correspondences)
            return None, count
    # no coloring possible or some edge has no correspondences
    if len(edges) == 0:
        return correspondences
    # still edges to fill in the correspondences for
    current_edge = edges[0]
    u,v = current_edge
    if len(correspondences[current_edge]) < min(colors[u], colors[v]) - 1:
        # current_edge is not going to be full
        remaining = edges[1:] + (current_edge,)
    else:
        # current edge will be full
        remaining = edges[1:]
    new_correspondences = dict(correspondences)
    if depth < branch_depth: # keep going, keep counting
        count = count+1 # for the current node
        for perm in add_correspondence_loop(correspondences[current_edge], colors[current_edge[0]]):
            if verbosity > 2:
                print("Let", current_edge, "have correspondence", perm)
            new_correspondences[current_edge] = perm
            test_bad_correspondence, count = bad_correspondence(tree, nontree, new_correspondences, remaining, colors, depth+1, count)
            if test_bad_correspondence is not None:
                return test_bad_correspondence, count
    elif depth == branch_depth and count % total_jobs == job_number: # keep going, don't count below here
        count = count+1 # for the current node
        for perm in add_correspondence_loop(correspondences[current_edge], colors[current_edge[0]]):
            if verbosity > 2:
                print("Let", current_edge, "have correspondence", perm)
            new_correspondences[current_edge] = perm
            test_bad_correspondence, _ = bad_correspondence(tree, nontree, new_correspondences, remaining, colors, depth+1, 0)
            if test_bad_correspondence is not None:
                return test_bad_correspondence, count
    elif depth > branch_depth: # keep going if you got here, but don't count
        for perm in add_correspondence_loop(correspondences[current_edge], colors[current_edge[0]]):
            if verbosity > 2:
                print("Let", current_edge, "have correspondence", perm)
            new_correspondences[current_edge] = perm
            test_bad_correspondence, _ = bad_correspondence(tree, nontree, new_correspondences, remaining, colors, depth+1, 0)
            if test_bad_correspondence is not None:
                return test_bad_correspondence, count
    # we get here if no bad_correspondence found or if branch_depth hit but we aren't continuing
    return None, count

def bad_correspondence_init(tree, nontree, colors):
    correspondences = dict()
    for e in nontree:
        correspondences[e] = ()
    return bad_correspondence(tree, nontree, correspondences, nontree, colors)

def all_colors_same(tree, nontree, num):
    res = dict()
    for i in vertices(tree, nontree):
        res[i] = num
    return res

# vertices is list of vertices
# edges is list of edges
# this function may need to wait for C version of this program
def spanning_trees_loop(vertices, edges):
    pass

# helper function
def gen_spanning_tree_recurse(graph, verts_in_tree, vert_to_expand, parent):
    nbs = graph[vert_to_expand]
    nontree = []
    subtrees = []
    added_nbs = []
    for j in nbs:
        if j in verts_in_tree: 
            if j is not parent: # adding this edge would make a loop, so it must a nontree edge
                nontree.append((vert_to_expand, j))
        else: # greedily adds all nbs to tree unless doing so would create a loop
            added_nbs.append(j)
            verts_in_tree.append(j)
    for j in nbs:
        if j is not parent and j in added_nbs: # we need to recurse on all the nbs which are now part of the tree
            subtree_j, nontree_j = gen_spanning_tree_recurse(graph, verts_in_tree, j, vert_to_expand)
            nontree += nontree_j
            subtrees.append(subtree_j)
    tree = (vert_to_expand, tuple(subtrees))
    return tree, nontree

# graph is dictionary with vertices as key, and values are list of neighbors
# output is a spanning tree and list of nontree edges suitable for plugging right into bad_correspondence_init
# function does not check for graph to be connected. If not so, tree and nontree will only cover one component
def gen_spanning_tree(graph):
    #global branch_depth
    verts_in_tree = [0] # always has 0 as root of tree
    tree, nontree = gen_spanning_tree_recurse(graph, verts_in_tree, 0, None)
    # nontree will contain each edge both ways, filter that out
    removes = []
    for pair in nontree:
        if pair in removes:
            continue
        else:
            u, v = pair
            removes.append((v,u))
    nontree_clean = tuple(e for e in nontree if e not in removes)
    #branch_depth = len(nontree_clean)
    return tree, nontree_clean

# following function has not been updated to handle spanning trees
#def Mycielski(graph):
    ## assumes vertices of graph are numbered 0 through n-1
    #n = len(vertices(graph))
    #res = list(graph)
    #for i, j in graph:
        #res.append((i, n+j))
        #res.append((i+n, j))
    #for k in range(n, n+n):
        #res.append((k, n+n))
    #return tuple(res)

#print(Mycielski(Mycielski(((0,1),))))

#G = Mycielski(Mycielski(Mycielski(((0,1),))))

#print(bad_correspondence_init(G, all_colors_same(G, 4)))

# example: K4 - e, or C4 + e, its the random graph!
#G_star_t = (1, ((2, ()), (3, ()), (4, ())))
#G_star_nt = ((2,4), (3,4))
#G_path_t = (1, ((3, ()), (2, ((4, ()),))))
#G_path_nt = ((1,4), (3,4))
#G_zig_t = (1, ((3, ()), (4, ((2, ()),))))
#G_zig_nt = ((1,2), (3,4))
#
#print("G is K4 - e or C4 + e.")
#print("Using star tree:")
#print(bad_correspondence_init(G_star_t, G_star_nt, all_colors_same(G_star_t, G_star_nt, 3)))
#print("Using path tree:")
#print(bad_correspondence_init(G_path_t, G_path_nt, all_colors_same(G_path_t, G_path_nt, 3)))
#print("Using zigzag tree:")
#print(bad_correspondence_init(G_zig_t, G_zig_nt, all_colors_same(G_zig_t, G_zig_nt, 3)))
#
#
## example K5 - e, which is still planar, so it should be 5 correspondence colorable, but is it 4? (Should be by Brooks' type result) Is it 3?
#H_star = [(1, ((2, ()), (3, ()), (4, ()), (5, ()))), ((2,3), (3,4), (4,5), (2,5), (3,5))]
#
#print("H is K5 - e")
#print("Using star tree:")
#print(bad_correspondence_init(H_star[0], H_star[1], all_colors_same(H_star[0], H_star[1], 4)))

verbosity = 1
n = 7

num_considered = 0
#branch_depth = 10 #n*(n-1)/2.0 - (n-1) - 5
job_number = sys.argv[1]
total_jobs = sys.argv[2]

for graph, graph6_str, choos, at_num, ub in graph6_file(f"choos_{n}.txt"):
    tree, nontree = gen_spanning_tree(graph)
    if verbosity > 0:
        print(f"Graph: {graph}")
    if verbosity > 1:
        print(f"Tree: {tree}, Nontree: {nontree}")
    at_num = int(at_num)
    k = int(choos)
    ub = int(ub)
    while k < ub: # when k == ub, we know it is correspondence coloring number
        branch_depth = math.ceil(math.log(6*total_jobs, k))
        thing = bad_correspondence_init(tree, nontree, all_colors_same(tree, nontree, k))
        if thing is None:
            break
        if verbosity > 0:
            print(f"Bad Correspondence for {k} colors: {thing}")
        k += 1
    if verbosity > 0:
        print(f"Correspondence coloring number {k}")
        print("")
    print(f"{graph6_str} {at_num} {k}")
