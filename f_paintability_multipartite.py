
from copy import copy, deepcopy
from functools import lru_cache
import itertools as it

def subsets(n):
    for k in range(n+1):
        for l in it.combinations(range(n),k):
            yield l

def lister(fvals):
    for t in it.product(*map(subsets, [len(x) for x in fvals])):
        # need to check for empty choice (not allowed)
        if all(len(l) == 0 for l in t):
            continue
        # otherwise it is fine
        yield t

def painter(fvals, lst):
    for i in range(len(fvals)):
        # can't pick this part if there are no vertices in it
        if len(lst[i]) == 0:
            #print(indent, "No vertices in this part to pick")
            continue
        new_fvals = [list(x) for x in fvals]
        # verts lister picked in this part are removed
        new_fvals[i] = [x for (j,x) in enumerate(fvals[i]) if j not in lst[i]]
        # verts lister picked in other parts have f values reduced by 1
        for j in range(len(fvals)):
            if j == i:
                continue
            for k in lst[j]:
                new_fvals[j][k] -= 1
        for j in range(len(fvals)):
            new_fvals[j].sort()
        tup_fvals = tuple([tuple(x) for x in new_fvals])
        yield tup_fvals

@lru_cache()
def f_paintable_multipartite(fvals):
    #print(indent, "fvals = ", fvals)
    # if there are no more vertices, then we can trivially color
    if len(fvals) == 0 or all(len(x) == 0 for x in fvals):
        #print(indent, "Trivial color")
        return True
    # if any value in fvals is 0, we cannot color
    if any(x == 0 for l in fvals for x in l):
        #print(indent, "Can't color")
        return False
    # otherwise lister picks a subset
    for lst in lister(fvals):
        #print(indent, "Lister picks: ", lst)
        val = False
        # painter picks an independent set (assume it is all of some part)
        for new_fvals in painter(fvals, lst):
            #print(indent, "Painter tries:", new_fvals)
            # recursively call function on new smaller graph
            # if it is paintable, then painter has a response to this choice of lister
            # no need to continue checking
            if f_paintable_multipartite(new_fvals):
                #print(indent, "Painter succeeds with:", new_fvals)
                val = True
                break
        # if painter has no good response, then lister wins
        if not val:
            #print(indent, "Lister wins this one")
            return False
    # if painter has a good response to all of lister's choices, painter wins
    #print(indent, "Painter wins this one")
    return True

def k_for_all_verts(k, parts):
    res = []
    for p in parts:
        res.append(tuple(it.repeat(k,p)))
    return tuple(res)

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


print(paintability((4,6),5))
