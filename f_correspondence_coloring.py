#!/usr/bin/python3

# Represent graph as list of edges, which are pairs of numbers (want to include number of verts?)

# Represent partial correspondence as list of lists (or tuples), where the ith element is what i is mapped to. Since we fill in the perms from 1 to n, this will work.

from functools import lru_cache
from itertools import count, filterfalse

num_considered = 0

# helper function
@lru_cache(maxsize=None)
def neighbors(v, graph):
    res = []
    for x,y in graph:
        if x == v:
            res.append((y, -1))
        elif y == v:
            res.append((x, 1))
    return res

# helper function
@lru_cache(maxsize=None)
def vertices(graph):
    res = []
    for x,y in graph:
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
    else: # orient == -1
        if c_uv in perm:
            return everything_but(perm.index(c_uv), k) # everything except what c_uv corresponds to
        else: # things that correspond to something
            return list(range(len(perm)))

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
    else: # orient == -1
        if c_v < len(perm): # c_v corresponds to something
            return perm[c_v] != c_u # so it can't be c_u
        else: # c_v does not correspond to anything 
            return c_u in perm # so c_u should correspond to something


# takes graph, partial correspondence, partial coloring, and list of uncolored vertex and attempts to extend coloring to the uncolored vertices
# graph: list of edges as pairs
# partial correspondence: dict of tuples, keys are edges, values are tuples giving partial correspondence for that edge
# partial coloring: dict of number, keys are vertices (also numbers), values are color of that vertex
# uncolored: list of uncolored vertices
# colors: dictionary giving number of colors available at each vertex
# returns extended coloring if exists, else returns None
# algorithm works by attempting to extend coloring to first element of list, then recursively calling itself on rest of list
def extend_coloring(graph, correspondences, coloring, uncolored, colors):
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
        for u, orient in neighbors(vertex, graph):
            if u in coloring:
                if orient == 1:
                    if not color_test(coloring[u], i, orient, correspondences[(u, vertex)]):
                        good_color = False
                        break
                else: # orient == -1
                    if not color_test(i, coloring[u], orient, correspondences[(vertex, u)]):
                        good_color = False
                        break
        if good_color:
            # i is good color for vertex
            new_coloring[vertex] = i
            # recursively color rest of graph
            extended_coloring = extend_coloring(graph, correspondences, new_coloring, remaining, colors)
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

def bad_correspondence(graph, correspondences, edges, colors):
    #print("Considering correspondence", correspondences)
    # shortcut if some edge of correspondences is ()
    if all(partial != () for partial in correspondences):
        test_coloring = extend_coloring(graph, correspondences, dict(), vertices(graph), colors)
        if test_coloring is not None:
            print("Coloring", test_coloring, "for correspondences", correspondences)
            return None
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
    for perm in add_correspondence_loop(correspondences[current_edge], colors):
        print("Let", current_edge, "have correspondence", perm)
        new_correspondences[current_edge] = perm
        test_bad_correspondence = bad_correspondence(graph, new_correspondences, remaining, colors)
        if test_bad_correspondence is not None:
            return test_bad_correspondence
    return None

def bad_correspondence_init(graph, colors):
    correspondences = dict()
    for e in graph:
        correspondences[e] = ()
    return bad_correspondence(graph, correspondences, graph, colors)

def all_colors_same(graph, num):
    res = dict()
    for i in vertices(graph):
        res[i] = num
    return res

def Mycielski(graph):
    # assumes vertices of graph are numbered 0 through n-1
    n = len(vertices(graph))
    res = list(graph)
    for i, j in graph:
        res.append((i, n+j))
        res.append((i+n, j))
    for k in range(n, n+n):
        res.append((k, n+n))
    return tuple(res)

#print(Mycielski(Mycielski(((0,1),))))

G = Mycielski(Mycielski(Mycielski(((0,1),))))


print(bad_correspondence_init(G, all_colors_same(G, 4)))

