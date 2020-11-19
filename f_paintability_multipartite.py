
from copy import copy, deepcopy
from functools import lru_cache
import itertools as it

def subsets(n):
    for k in range(n+1):
        for l in it.combinations(range(n),k):
            yield l

def upto(n):
    for i in range(n+1):
        yield i

def leqList(l):
    for k in it.product(*map(upto, l)):
        yield k

def lister(fvals):
    for t in it.product(*map(leqList, fvals)):
        # need to check for empty choice (not allowed)
        if all(sum(l) == 0 for l in t):
            continue
        # otherwise it is fine
        yield t

def decNotPicked(k, l):
    return tuple(k[i] - l[i] + l[i+1] if i < len(k)-1 else k[i] - l[i] for i in range(len(k)))

def decPicked(k, l):
    return tuple(k[i] - l[i] for i in range(len(k)))

def painter(fvals, lst):
    for i in range(len(fvals)):
        # can't pick this part if there are no vertices in it
        if sum(lst[i]) == 0:
            #print(indent, "No vertices in this part to pick")
            continue
        yield tuple(decNotPicked(fvals[j], lst[j]) if j != i else decPicked(fvals[j], lst[j]) for j in range(len(fvals)))

@lru_cache()
def f_paintable_multipartite(fvals, indent='', verbose=False):
    if verbose: print(indent, "fvals = ", fvals)
    # if a vertex has higher f value than its degree, we can delete it
    # (we can paint by only picking it if none of its neighbors are picked,
    # since f value > degree, this is guarenteed to happen. So it can be colored,
    # and we can ignore it.)
    n = sum(sum(l) for l in fvals) # number of vertices in graph
    while any(len(l) > n-sum(l)+1 for l in fvals): # recursively delete vertices until all f values at most degrees
        if verbose: print(indent, "Removing vertices with f value greater than degree")
        fvals = tuple(l[:n-sum(l)+1] for l in fvals) # +1 to allow f value to equal degree
        n = sum(sum(l) for l in fvals) # update number of vertices in graph
        if verbose: print(indent, "New fvals =", fvals)
    # if there are no more vertices, then we can trivially color
    if len(fvals) == 0 or all(sum(x) == 0 for x in fvals):
        if verbose: print(indent, "Trivial color")
        return True
    # if any vertices have an f value of 0, we cannot color
    if any(l[0] > 0 for l in fvals):
        if verbose: print(indent, "Can't color")
        return False
    # otherwise lister picks a subset
    for lst in lister(fvals):
        if verbose: print(indent, "Lister picks: ", lst)
        val = False
        # painter picks an independent set (assume it is all of some part)
        for new_fvals in painter(fvals, lst):
            if verbose: print(indent, "Painter tries:", new_fvals)
            # recursively call function on new smaller graph
            # if it is paintable, then painter has a response to this choice of lister
            # no need to continue checking
            if f_paintable_multipartite(new_fvals, indent=indent+"\t", verbose=verbose):
                if verbose: print(indent, "Painter succeeds with:", new_fvals)
                val = True
                break
        # if painter has no good response, then lister wins
        if not val:
            if verbose: print(indent, "Lister wins this one")
            return False
    # if painter has a good response to all of lister's choices, painter wins
    if verbose: print(indent, "Painter wins this one")
    return True

def k_for_all_verts(k, parts):
    return tuple(map(lambda x: (0,)*k + (x,), parts))

def paintability(parts, kmax):
    for k in range(2, kmax):
        if f_paintable_multipartite(k_for_all_verts(k, parts)):
            return k
    return -1

def test_k(k, imin, imax, jmin, jmax):
    for i in range(imin,imax):
        for j in range(max(i, jmin),jmax):
            print("Complete", i, ",", j)
            parts = (i,j)
            print(f_paintable_multipartite(k_for_all_verts(k, parts)))


test_k(3, 2, 8, 2, 8)
