#------------------------------------------------------
# Contents
#------------------------------------------------------

#    rev_powerset
#    SubgraphsDealer
#    CNS
#    IndSets
#    restriction_list_single_check
#    set2tup
#    reduce
#    timestring
#    subsetstring
#    fChoosable
#    restrictedListAssignments
#    assignmentCheck




#------------------------------------------------------
# Libraries, Comments, Global Variables
#------------------------------------------------------

# These scripts were written with the intent of using Python 2.7 before the switch to Python 3.  There should not be many differences to implement our code in Python 3 once this change is made.

import copy as cpy
from sage.all import *
from sage.graphs.graph import Graph
from sage.graphs.generators.basic import CompleteGraph
import time
from sage.rings.polynomial.polynomial_ring_constructor import *
from sage.rings.finite_rings.finite_field_constructor import *
ZZ = sage.rings.integer_ring.IntegerRing()
from sage.graphs.independent_sets import IndependentSets
from itertools import combinations
from numpy import prod





#------------------------------------------------------
# The Routines
#------------------------------------------------------



def rev_powerset(n):  # reverse ordering of colex on sets.
    if n > 0:
        for subset in rev_powerset(n-1):
            subset.append(n-1)
            yield subset
        for subset in rev_powerset(n-1):
            yield subset
    else:
        yield []
    return






class SubgraphsDealer:
    
    def __init__(self,subgraphs_list,LAscheme="E"):
        self.subgraphs_list = cpy.copy(subgraphs_list)
        self.LAscheme = LAscheme
        
        # subgraph is the current connected subgraph (set of vertices).
        self.subgraph = None
        
        self.required_vertex = None
    
    
    
    def copy(self):
        clone = SubgraphsDealer(cpy.copy(self.subgraphs_list),self.LAscheme)
        clone.subgraph = cpy.copy(self.subgraph)
        return clone
    
    
    
    def prune_subgraphs_list_with_E(self,E):
        i = 0
        while i < len(self.subgraphs_list):
            A = self.subgraphs_list[i]
            # if A not subset E, then toss A from the list
            if not A <= set(E):
                del self.subgraphs_list[i]
                continue
            else:
                i += 1
    
    
    # call if small condition met
    def prune_subgraphs_list_with_CNS(self,G_copy, ordering, f, CA):
        i = 0
        while i < len(self.subgraphs_list):
            A = self.subgraphs_list[i]
            # if A can be be shown to be unnecessary via CNS, then toss A from the list
            if restriction_list_single_check(G_copy, ordering, f, CA, A):
                del self.subgraphs_list[i]
                continue
            else:
                i += 1
    
    
    
    
    # E should be nonempty to call this
    def choose_required_vertex(self,L,E):
        
        if self.LAscheme == "E":
            self.required_vertex = E[0]
        
        elif self.LAscheme == "L":
            m = min([(L[y],y) for y in E])
            self.required_vertex = m[1]
    
    
    
    # before calling this, call choose root.
    # before calling this for the first time, call prune by E (and possibly by CNS)
    def next_with_required_vertex(self):
        # go through list to find A containing required_vertex
        i = 0
        while i < len(self.subgraphs_list):
            # if required_vertex not in A, then keep A in the list
            A = self.subgraphs_list[i]
            if not self.required_vertex in A:
                i += 1
                continue
            # if A works, then set .subgraph to A, remove A, and return True
            else:
                self.subgraph = A
                del self.subgraphs_list[i]
                return True
        # if no such A, then return False (v can't be filled)
        else:
            return False
    
    
    # before calling this, call choose root
    def next_with_required_vertex_not_containing_set(self,D):
        while self.next_with_required_vertex():
            if D <= self.subgraph:
                continue  # try again -- we can toss that subgraph
            else:
                return True
        else:
            return False  # no subgraph remaining
    
    
    




















def CNS(G, f, verbosity=0):
    
    if verbosity > 7:
        print_flag = True
    else:
        print_flag = False
    
    # INITIALIZE
    
    # n is the number of vertices.
    n = G.order()
    if n != len(f):
        print "Error:  Length of f does not match order of G."
        return

    if [x for x in f if x<1]:
        if print_flag:  print "Error: f-vector contains non-positive entries."
        return False
    
    # gapleft is the remaining number of steps below the sum of (f-1)-values for our monomial.
    # gapleft begins as the disparity and should drop down to 0.
    # In general, it is given by disparity - sum([f[x]-1-deg[x] for x in range(v)]).
    gapleft = sum(f) - len(f) - G.size()
    if gapleft < 0:
        if print_flag:  print "    >> Disparity:  %d.  Cannot apply Nullstellensatz."%(gapleft)
        return False
    elif print_flag:  print "    >> Disparity:  %d.  Proceeding."%(gapleft)
    
    # Relabel graph so that expansion steps have smaller memory requirements.
    perminv = []
    ma = G.order()
    while len(perminv) < ma:
        degrees = [len(set(G.neighbors(x))-set(perminv)) for x in G.vertices()]
        for j in perminv:
            degrees[j] = ma
        mi = min(degrees)
        i = degrees.index(mi)
        perminv.append(i)
    
    if print_flag:  print "    >> Relabeling via",perminv
    # Want perminv[i] to be relabeled as i.
    perm = [perminv.index(x) for x in range(len(perminv))]
    G.relabel(perm)
    f = [f[x] for x in perminv]
    
    # Polynomial stuff.
    stringlist = ['x'+str(x) for x in range(n)]
    P=PolynomialRing(ZZ, n, stringlist)
    varlist = P.objgens()[1]
    
    # v is the vertex/variable we are inspecting.
    v = 0
    
    # deg, mins, and coeff all act as stacks.
    deg = []    # deg is a list holding the degrees we're currently using.
    mins = []   # mins is a list holding the lower bounds for the degrees we can use in deg.
    coeff = []  # coeff is a list holding the coefficients we've found for varlist[i]**deg[i].
    # Entries in deg start high, but go down by 1 to the corresponding mins value as we backtrack.
    
    # The coeff stack is always one longer than deg and mins; the first entry should be the polynomial 1.
    coeff.append(varlist[0]**0)
    
    
    # START LOOP
    
    while True:
        
        if coeff[-1] == 0 or gapleft < 0:  # If we have already encountered an issue with expansion.
            
            # BACKTRACK
            # If the current monomial expansion has already resulted in a coefficient of 0,
            # then the coefficient is 0 in the entire polynomial and we need to backtrack.
            # Lower the previous variable's degree to see if we get higher degrees afterward.
            
            while True:
                
                del coeff[-1]
                v -= 1
                
                # If the current degree is more than the minimum, we can break out of the loop to drop down now.
                if deg[-1] > mins[-1]:
                    break
                # Otherwise, march backwards.  We'll keep going until we can drop without going below mins.
                
                gapleft += f[v]-1-deg[-1]
                del deg[-1]
                del mins[-1]
                
                # If we've marched back through everything, then we've exhausted our options.
                if v < 1:
                    if print_flag:  print "    >> The CNS is inconclusive; no relevant monomials have nonzero coefficient."
                    G.relabel(perminv)
                    f = [f[x] for x in perm]
                    return False
            
            # Drop down the degree (and note this in diff).
            deg[-1] -= 1
            gapleft -= 1
            # Now that we've backtracked, we try moving forward again by finding the next coefficient piece.
        
        
        else:  # If there is hope for expanding the current monomial.
            
            # FORWARD
            
            # If we have built the entire monomial, then stop.
            if v > n-1:
                if print_flag:  print "    >> The graph is choosable!  \n    >>Certificate:  The monomial with power list "+str(deg)+" has coefficient "+str(coeff[-1])+"."
                G.relabel(perminv)
                f = [f[x] for x in perm]
                return True
            
            # Otherwise, append next values to deg and mins.
            
            # f[v]-1-gapleft is a lower bound for deg[v].  (Any power less than this number 
            # will cause the monomial to already have too low of total degree, even if the 
            # remaining powers are all the f-1 values.)
            mins.append(max(0,f[v]-1-gapleft))
            
            a = coeff[-1].degree(varlist[v])  # a is the degree of xv in the residual coefficient.
            b = len([x for x in G.neighbors(v) if v < x])  # b is the number of forward neighbors of v.
            
            # The monomial we have built so far cannot have a power of xv higher than a+b.
            if a+b < f[v]-1:
                deg.append(a+b)
                gapleft -= f[v]-1-a-b
            else:
                deg.append(f[v]-1)
            # Now find the next coefficient piece and try to keep going.
        
        
        # FIND THE NEXT COEFFICIENT PIECE
        
        # Expand remaining polynomial in the vth variable (times the previous coefficient, which might contain 
        # the vth variable) to find coefficient of varlist[:v+1]**deg[:v+1].
        
        prod = varlist[0]**0
        for factor in [varlist[v]-varlist[u] for u in G.neighbors(v) if v < u]:
            prod *= factor
        
        newco = 0
        for x in range(0,deg[-1]+1):
            newco += coeff[-1].coefficient({varlist[v]:x}) * prod.coefficient({varlist[v]:(deg[-1]-x)})
        coeff.append(newco)
        
        
        # Increment and go back to the top of the loop!
        v += 1
        






# For using the .next() method.
def IndSets(graph,maximal=True,complement=False):
    for I in IndependentSets(graph,maximal=maximal,complement=complement):
        yield I









# Checks whether the subset A is CNS-reducible when adding it to the current list assignment CA.  Returns True if the subset is reducible for the given CA; returns False otherwise (which means it must stay in the restriction list).
def restriction_list_single_check(G_copy, ordering, f, CA, A):
    # we test colorings from CA.append(A) by using backtracking with a stack of color classes
    ISgens = []      # elements of ISgens will be IndependentSets iterators
    coloring = []    # elements of coloring will be pairs (maximal independent set, union of all first entries up to this pair)
    c = 0            # at the top of the following loop, c is the color we are about to assign a color class (an IndependentSets item)
    while True:
        if c==len(CA):
            ordering_template = list(ordering)
            newf_template = list(f)
            for C in CA:
                for v in C:
                    newf_template[v] -= 1
            subG_template = Graph(G_copy)
            S = A    # S is going to be A, but without vertices which have already been colored
            if c>0:
                S = S - coloring[c-1][1]
                for v in coloring[c-1][1]:
                    del newf_template[ordering_template.index(v)]
                    del ordering_template[ordering_template.index(v)]
                    subG_template.delete_vertex(v)
            H = Graph(G_copy)
            for v in H.vertices():
                if v not in S:
                    H.delete_vertex(v)
            disparity_temp = sum(newf_template) - len(newf_template) - subG_template.size()
            # even if the disparity is negative, coloring certain independent sets in certain classes might still make the CNS applicable to the remaining graph.
            for I in IndSets(H):  # check if the coloring given by 'coloring' and I shows via CNS that A is reducible
                if sum([newf_template[ordering_template.index(v)]-1 for v in I]) + len(S)-len(I) - sum([subG_template.degree(v) for v in I]) <= disparity_temp:
                    newf = list(newf_template)
                    newordering = list(ordering_template)
                    subG = Graph(subG_template)
                    for v in I:
                        del newf[newordering.index(v)]
                        del newordering[newordering.index(v)]
                        subG.delete_vertex(v)
                    for v in S-set(I):
                        newf[newordering.index(v)] -= 1
                    subG.relabel()
                    newf = [newf[newordering.index(x)] for x in range(G_copy.order()) if x in newordering]

                    subG,newf,same,low = reduce(subG,newf)
                    
                    if low:  # In this case, the remaining subG,newf are NOT f-choosable.  Try another I.
                        continue

                    if len(newf) < 1 or CNS(subG,newf):
                        return (set(I),A)

            # otherwise, backtrack as far as necessary and take one step forward
            while c>0:
                try:
                    del coloring[c-1]
                    new = set(ISgens[c-1].next())  # this throws the exception if the generator is exhausted
                    if c==1:
                        coloring.append((new,new))
                    else:
                        coloring.append((new,new|coloring[c-2][1]))
                    break
                except StopIteration:
                    del ISgens[c-1]
                    c -= 1
            if c==0:
                break
            else:
                continue

        # if c < len(CA), we come here and add a fresh generator to our stack
        if c==0:
            S = CA[0]
        else:
            S = CA[c] - coloring[c-1][1]
        H = Graph(G_copy)
        for v in H.vertices():
            if not v in S:
                H.delete_vertex(v)
        ISgens.append(IndSets(H,maximal=True))  # Using maximal works since they are maximal with respect to the *remaining* graph
        new = set(ISgens[c].next())
        if c==0:
            coloring.append((new,new))
        else:
            coloring.append((new,new|coloring[c-1][1]))
        c+=1

    return False
























# To handle the unorderedness of sets.  Takes a set S of nonnegative integers and returns a tuple of the same elements, ordered by <.  Could easily be adjusted to respect ordering by passing ordering in.
def set2tup(S):
    return tuple(x for x in range(max(S)+1) if x in S)






# Takes G,f as input and returns the pruned version, having removing any vertices with f-values higher than degree (iteratively) and having removed precolored vertices (f-values of 1) and adjusted the neighborhoods accordingly (also iteratively).
def reduce(G_copy,f_copy):
    G = Graph(G_copy)
    f = f_copy[:]
    redundance = [x for x in G.vertices() if f[x] > G.degree(x)]
    while redundance:
        for v in redundance:
            G.delete_vertex(v)
        G.relabel()
        for v in reversed(redundance):
            del f[v]
        redundance = [x for x in G.vertices() if f[x] > G.degree(x)]
    ones = [x for x in G.vertices() if f[x]==1]
    while ones:
        v = ones[0]
        for u in G.neighbors(v):
            f[u] -= 1
        del f[v]
        G.delete_vertex(v)
        G.relabel()
        ones = [x for x in G.vertices() if f[x]==1]
    same = (len(f_copy)==len(f)) # same indicates whether the returned configuration is the same as the input configuration
    low = (len([x for x in f if x<1])>0) # low indicates whether there are nonpositive f-values in the output configuration
    return G,f,same,low









def timestring(time,setting="m"):
    t = int(time)
    d = str(int((10*(time-t))%10))
    if setting=="s":
        string = str(t)+"."+d+"s"
    else:
        s = str(t%60)
        m = int(t/60)
        if setting=="m":
            string = str(m)+"m "+s+"."+d+"s"
        else:
            h = int(m/60)
            m %= 60
            if setting=="h":
                string = str(h)+"h "+str(m)+"m "+s+"."+d+"s"
            else:
                d = int(h/24)
                h %= 24
                string = str(d)+"d "+str(h)+"h "+str(m)+"m "+s+"."+d+"s"
    return string







def subsetstring(A,n):
    string = " "*8+"["
    for i in range(n):
        string += " "*3
        if i in A:
            string += str(i)
        else:
            string += " "*len(str(i))
    string += "   ]"
    return string
















def fChoosable(G_copy, f_copy, inprocess=3, print_mod=1000, LAscheme="E", rev=False, verbosity=0, search_all=0):
    
    if not LAscheme in ["L","E"]:
        print 'Error:  LAscheme must be one of "L" or "E".'
        return

    if verbosity>0:
        start = time.clock()
        num_pruned = 0

    if len(f_copy) != G_copy.order():
        print "Error: f-vector length does not match graph order."
        return

    if [x for x in f_copy if x<1]:
        print "Error: f-vector contains non-positive entries.  (If correct f, graph is not f-choosable.)"
        return

    G,f,same,low = reduce(G_copy,f_copy)
    
    if not same:
        if verbosity>0:  print "The input configuration was reduced for having f-values either higher than the degree or equal to 1."
        if not f:
            if verbosity>0:  print "In fact, it reduced entirely; the graph is f-choosable via a greedy argument."
            return (True,'greedy')
        if verbosity>0:  print "f:",f
    if low:
        if verbosity>0:  print "Reduced f-vector contains non-positive entries.  Graph is not f-choosable."
        return (False,'precoloring')
    
    if verbosity>0:  print "Testing:  Is the graph (%d vertices, %d edges) "%(G.order(),G.size())+str(f)+"-choosable? \n\nFirst, we consider whether the Combinatorial Nullstellensatz applies."
    if verbosity>0:  boo = CNS(G,f,verbosity=10)
    else:  boo = CNS(G,f,verbosity=0)
    if verbosity>0:  print "    (Time so far:  "+timestring(time.clock()-start)+")"
    if boo:
        return (True,'CNS')

    if verbosity>0:  start2 = time.clock()

    if verbosity>0:  print "\nNow we move on to our exhaustive search.  \nTo start up our process, we find and remove connected subgraphs which cannot be part of a bad list assignment, and then check the multiplicities of those that remain."

    it = restrictedListAssignments(G, f, inprocess_depth=inprocess, LAscheme=LAscheme, rev=rev, verbosity=verbosity)
    
    bad_list_assignments = []

    if verbosity>0:
        nodes_count = 0
        full_count = 0
        bottom_count = 0
        node_code = [0,]
        bad_LA_count = 0

    try:
        x = it.next()           # The first call to a Python iterator must be .next().
        # ^^Recall: x is a tuple (colClasses,E_stack,generator_stack,back,num_pruned).
                
        if verbosity>0:
            print "    >> Start-up complete.  Total time so far:  "+timestring(time.clock()-start)
        
        
            print "\nNow, begin building and checking list assignments with the in-process pruning parameter set to %d. \nInformation will print every %d nodes we visit in the search space."%(inprocess,print_mod)
            print "( a | b ) means that colorability class is the a-th of b subgraphs on that level of the stack."
            start3 = time.clock()
            lap = start3

            nodes_count += 1
            node_code[0] += 1

        # At the top of this loop, we have obtained an x from the listAssignment iterator.
        while True:

            if not x[1][-1] and verbosity>0:
                full_count += 1

            c = assignmentCheck(G,x[0])

            # If the list assignment is not colorable and full, then report it and stop.
            if not c and not x[1][-1]:
                if verbosity>0:
                    bad_LA_count += 1
                    print " "*10+"*"*10+" BAD LIST ASSIGNMENT # %d:"%(bad_LA_count)+"*"*10
                    for i in range(len(x[0])):
                        print " "*10+subsetstring(x[0][i],G.order()),"(",node_code[i],")"# |",x[3][i],")"
                    print
                if search_all>0:
                    # Then switch c to be nonempty so that the iterator backtracks.
                    c = [0,]
                    if search_all>1:
                        bad_list_assignments.append(cpy.copy(x[0]))
                else:
                    bad_list_assignments.append(cpy.copy(x[0]))
                    break


            if verbosity>0 and x[3]:
                bottom_count += 1
            
            # Tell the listAssignment iterator to go to the next relevant list assignment (either moving forward or backtracking).
            x = it.send(c)
            
            if verbosity>0:
                if len(x[0]) > len(node_code):
                    node_code.append(1)
                else:
                    node_code[len(x[0])-1] += 1
                    for i in range(len(x[0]),len(node_code)):
                        del node_code[-1]

                num_pruned += x[4]

                nodes_count += 1
                if nodes_count%print_mod == 0:
                    print "-----> ",nodes_count,"partial list assignments inspected so far.  # full LAs:",full_count,".  # bottom LAs:",bottom_count,".  # bad LAs:",bad_LA_count,".  # inprocess prunes:",num_pruned,"."
                    print "^^^^^>  Time taken on this batch of PLAs:",timestring(time.clock()-lap),"    Time taken on the whole (subgraph) job so far:",timestring(time.clock()-start)
                    print "^^^^^>  Current partial list assignment:"
                    print subsetstring(x[0][0],G.order()), "(", node_code[0], "|", len(x[2][0].subgraphs_list), ")"
                    for i in range(1,len(x[0])):
                        print subsetstring(x[0][i],G.order()), "(", node_code[i], "|", len(x[2][i].subgraphs_list), ")"
                    print "\n"*3
                    lap = time.clock()

    except StopIteration:
        if verbosity>0:  print "\nNo bad list assignments!\n"

    if verbosity>0:
        if nodes_count < 1:  # If CNS inconclusive but a vertex is not contained by any remaining subgraphs, then the normal start3 gets skipped.
            start3 = time.clock()
        end = time.clock()
        print "Finished.  %d PLAs visited which weren't pruned via CNS methods."%(nodes_count)
        print "# full LAs:  %d.  # bottom LAs:  %d.  # bad LAs:  %d.  # inprocess prunes:  %d"%(full_count,bottom_count,bad_LA_count,num_pruned)
        print "Time taken on initial CNS:       "+timestring(start2-start)
        print "Time taken on subgraph pruning:  "+timestring(start3-start2)
        print "Time taken on exhaustive search: "+timestring(end-start3)
        print "Total time taken on entire job:  "+timestring(end-start)
    if not bad_list_assignments:
        return (True,'brute')
    elif search_all>1:  # stored all bad list assignments
        return (False,'brute',bad_list_assignments)
    elif search_all<1:  # found one bad list assignment and then stopped
        return (False,'brute',bad_list_assignments[0])
    else:  # printed all bad list assignments
        return (False,'brute')





# inprocess_depth controls how many partial list assignments are checked for the reducibility when extending during the normal routine.
def restrictedListAssignments(G,f,inprocess_depth=3,LAscheme="E",rev=False,verbosity=0):

    ordering = range(G.order())
    reducers=[]

    # Initialize the things!
    generator_stack = []
    colClasses = []  # colClasses[i] is {vertices with color i in their list}.  We yield this.
    L = [0]*G.order()  # Keeps track of how many colors each vertex is currently assigned.
    D = None
    E_stack = []
    E_stack.append([x for x in ordering if L[x]<f[x]])  # Eligible vertices for colorability classes.
    col = 0  # Keeps track of next colorability class('s color).
    ticket = None

    reducers = [set(x) for x in reducers]
    
    
    
    CA = []
    
    if verbosity>1:
        lap = time.clock()
        print "    >> Now vetting all connected subgraphs."

    initial_restriction_list=[]
    
    reducibleList = []
    
    if verbosity>1:
        tenpercent = int(2**G.order()/100)+1
        counter = 0
        lap = time.clock()
        begin = lap
    
    it = rev_powerset(G.order())
    A = set(it.next())    
    while A:  # set up this way to avoid the empty set.
    
        if not G.subgraph(A).is_connected():
            A = set(it.next())
            continue
    
        if verbosity>1:
            if counter %tenpercent==tenpercent-1:
                print "        %d vetted.  Time on batch:  "%(counter)+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1

        red = False

        # if A is reducible for CA, continue/break to move on to next subset

        # first see if A is reducible by previous observations
        for x in reducibleList:
            if x[0] <= A and A <= x[1]:
                red = True
                break
        if red:
            A = set(it.next())
            continue

        # now see if A can be shown reducible by CNS
        red = restriction_list_single_check(G, ordering, f, CA, A)  # red is either False or a tuple (I,A) where I is the independent set that worked with the CNS to show that A is "reducible" with CA

        if not red:
            initial_restriction_list.append(A)
        else:
            reducibleList.append(red)
        
        A = set(it.next())

    if verbosity>1:
        print "        Time vetting all:  "+timestring(time.clock()-begin)
        print "    >>",len(initial_restriction_list),"remaining subgraphs can potentially be used to build list assignments."
        
        if len(initial_restriction_list) == 0:
            print "    >> So we can't build any bad list assignments!"
            return

        tenpercent = int(len(initial_restriction_list)/100)+1
        lap = time.clock()
        counter = 1
        begin = lap
        
    num_pruned = 0  # need to define num_pruned regardless of verbosity since it is passed in ticket.

    max_mult_init = {}
    for A in initial_restriction_list:
        
        CA = [A,]

        cap = min([f[x] for x in A])
        while len(CA) < cap:

            if restriction_list_single_check(G, ordering, f, CA, A):
                break

            CA.append(A)
        max_mult_init[set2tup(A)] = len(CA)
        
        if verbosity>1:
            if counter%tenpercent==0:
                print "       ",counter,"multiplicities checked.  Time on batch:  "+timestring(time.clock()-lap)
                lap = time.clock()
            counter += 1


    if verbosity>1:
        print "        Time for all multiplicity checks:  "+timestring(time.clock()-begin)
        print "    >> Multiplicities breakdown:"
        for i in range(1,max(max_mult_init.values())+1):
            print "        Number of subgraphs with multiplicity %d:  %d"%(i,len([x for x in max_mult_init.values() if x==i]))
        print "    >> Vertex containment breakdown for these leftover subgraphs:"
        for v in range(G.order()):
            print "        Number containing vertex %d: "%(v),len([x for x in initial_restriction_list if v in x])

    subgraph_pool = {frozenset(x) for x in initial_restriction_list}
    perminv = []
    ma = sum(max_mult_init.values())+1
    while len(perminv) < G.order():
        containments = [sum([max_mult_init[set2tup(x)] for x in subgraph_pool if v in x]) for v in G.vertices()]
        for j in perminv:
            containments[j] = ma
        mi = min(containments)
        i = containments.index(mi)
        perminv.append(i)
        batch = [x for x in subgraph_pool if i in x]
        subgraph_pool -= set(batch)
    
    if verbosity>1:  print "Relabeling everything using perminv:",perminv
    # Want perminv[i] to be relabeled as , and want j to be relabeled as perm[j].
    perm = [perminv.index(x) for x in range(len(perminv))]
    G.relabel(perm)
    f = [f[x] for x in perminv]
    relabeled_subgraph_pool = {frozenset([perm[i] for i in x]) for x in initial_restriction_list}

    max_mult = {}
    current_mult = {}
    for A in relabeled_subgraph_pool:  # A is new, but max_mult_init is for old labelings -- go back using perminv.
        max_mult[set2tup(A)] = max_mult_init[set2tup({perminv[a] for a in A})]
        current_mult[set2tup(A)] = 0

    subgraphs_list = []
    mult = 1
    while relabeled_subgraph_pool:
        subpool = {x for x in relabeled_subgraph_pool if max_mult[set2tup(x)]==mult}
        relabeled_subgraph_pool -= subpool
        size = G.order()
        while subpool:
            batch = [x for x in subpool if len(x)==size]
            subpool -= set(batch)
            for x in batch[::-1]:
                subgraphs_list.append(set(x))
            size -= 1
        mult += 1
    if rev:
        subgraphs_list = subgraphs_list[::-1]
        
    
    it = SubgraphsDealer(subgraphs_list,LAscheme=LAscheme)
    
    generator_stack.append(it)                         # Stack of subgraph iterators; cth element generates the c-colorability class.

    it.choose_required_vertex(L,E_stack[0])  # it.required_vertex = ordering[0]
    if it.next_with_required_vertex():
        A = generator_stack[col].subgraph  # generator_stack[col] == it
    
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])  # for x in E
    else:
        return

    back = False

    # Begin generation!
    while generator_stack:
        
        # Yield the current assignment; receive information from caller.
        ticket = yield (colClasses,E_stack,generator_stack,back,num_pruned)
        back = False
        
        # Case 1:  Move forward.  (Build the current partial list assignment up one more level.)
        # Since it is a new color, we start/copy a new iterator and add it to the stack (rather than call one we already have).
        # (There is no use for D in this case.)
        if not ticket:  # Case:  The partial list assignment was not colorable.
            # Add the new colorability class iterator to the iterator stack, then call the first subgraph, add it to the stack, and update.
            generator_stack.append(generator_stack[col].copy())  # Add new subgraph iterator (for the color col+1) to the stack.
            col+=1                                               # Increment col -- now the class we are generating is for the color col.

            # Because we've just made a copy of the top generator, we need to check the current colorability class for repetition:
            A = generator_stack[col].subgraph
            if current_mult[set2tup(A)] < max_mult[set2tup(A)] and A <= set(E_stack[-1]):
                if len(colClasses) < inprocess_depth:
                    if not restriction_list_single_check(G, ordering, f, colClasses, A):
                        colClasses.append(A)
                        current_mult[set2tup(A)] += 1
                        for v in colClasses[col]:                            # Update L.
                            L[v]+=1
                        E_stack.append([x for x in ordering if L[x]<f[x]])         # Update E.
                        continue
                else:# If we have more than a few colors and we know A fits, let's just try it!
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:                            # Update L.
                        L[v]+=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])         # Update E.
                    continue
            # If the current subgraph doesn't fit into E, then we check for the next eligible subgraph after updating restriction lists.
            # We don't update the restriction lists when backtracking since they've already been updated here first.
            
            
            generator_stack[col].choose_required_vertex(L,E_stack[-1])
            generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])

            # (Update the restrictionList from the current colClasses if we are in the early stages.)
            if len(colClasses) < inprocess_depth:
                generator_stack[col].prune_subgraphs_list_with_CNS(G, ordering, f, colClasses)
            
            if generator_stack[col].next_with_required_vertex():
                A = generator_stack[col].subgraph
                colClasses.append(A)
                current_mult[set2tup(A)] += 1
                for v in colClasses[col]:                            # Update L.
                    L[v]+=1
                E_stack.append([x for x in ordering if L[x]<f[x]])         # Update E.
                continue
            # If there is no such subgraph, then we are as full as possible and need to backtrack.
            back = True

        # Case 2a:  Backtrack because we tried to fill the current list assignment, but there were no more valid subsets.
        # Skip over the following else portion down to the main backtracking.


        # Case 2b:  Backtrack because the current list assignment colClasses is colorable (on the whole graph).
        # Start by moving the top color forward with respect to not containing the current colorability class.
        else:
            back = True
            # Delete the top colorability class (but not the generator for the class quite yet).
            # print col, colClasses[col], set2tup(colClasses[col]), colClasses
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            del E_stack[-1]

            # If the provided coloring (or lack of a provided coloring) is not on all the vertices, then we can't skip ahead in any clever way.
            if len(ticket) < G.order():
                generator_stack[col].choose_required_vertex(L,E_stack[-1])
                generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])
                if generator_stack[col].next_with_required_vertex():
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])
                    continue
                # If not, then we need to peel off more generators/layers.

            # But if the provided coloring is on all vertices, then we can use the top color class to potentially skip past some colorability class candidates.
            else:
                D = [x for x in colClasses[-1] if ticket[x]==col]  # D is the set of vertices in the coloring's highest color class.

                # First, see if the generator can give us something new:
                generator_stack[col].choose_required_vertex(L,E_stack[-1])
                generator_stack[col].prune_subgraphs_list_with_E(E_stack[-1])
                if generator_stack[col].next_with_required_vertex_not_containing_set(D):
                    A = generator_stack[col].subgraph
                    colClasses.append(A)
                    current_mult[set2tup(A)] += 1
                    for v in colClasses[col]:
                        L[v] +=1
                    E_stack.append([x for x in ordering if L[x]<f[x]])
                    continue
                # If not, then we need to peel off more generators/layers.  Note: D no longer serves any purpose past here.


        # At the moment, generator_stack[col] is an iterator that has no valid subgraphs left (either because it was a false start in Case 2a or because we backtracked through it in Case 2b).  But colClasses[col] is undefined (either because it didn't get there in Case 2a or because we deleted it in Case 2b).
        # Case 2:  Backtrack.  Step back as many colors as necessary and then obtain a next subgraph -- add that to colClasses.

        # At this point, we are guaranteed to take at least one step back through our stacks.
        del generator_stack[col]
        col -= 1
        if col<0:
            return
        for v in colClasses[col]:
            L[v] -=1
        current_mult[set2tup(colClasses[col])] -= 1
        del colClasses[col]
        del E_stack[-1]

        # Now step back as many times as necessary.
        generator_stack[col].choose_required_vertex(L,E_stack[-1])
        while not generator_stack[col].next_with_required_vertex():

            del generator_stack[col]
            col -= 1
            if col<0:  # If generator_stack is empty, then we have backtracked entirely.  A StopIteration error will be thrown.
                return
            # Prepare for bumping the now-top generator on the stack.
            for v in colClasses[col]:
                L[v] -=1
            current_mult[set2tup(colClasses[col])] -= 1
            del colClasses[col]
            del E_stack[-1]
            generator_stack[col].choose_required_vertex(L,E_stack[-1])

        # Now we have successfully bumped the generator to another subgraph.  Push the subgraph onto the colClass stack.
        A = generator_stack[col].subgraph
        colClasses.append(A)
        current_mult[set2tup(A)] += 1
        for v in colClasses[col]:
            L[v]+=1
        E_stack.append([x for x in ordering if L[x]<f[x]])
        continue















def assignmentCheck(G,CAtocopy):
    # Convert the colorability classes assignment to a list assignment.
    CA = []
    for a in CAtocopy:
        CA.append(cpy.deepcopy(a))
    n = G.order()
    LA = []                     # LA is the list assignment, where L[v] is the set of colors assignable to vertex v.
    for j in range(n):          # It will change according to our partial colorings during the algorithm.
        LA.append(set())
    for col in range(len(CA)):
        for v in CA[col]:
            LA[v].add(col)

    # Initialize!
    i = 0          # Which vertex we're coloring at the moment.
    c = [-1]*n     # Our coloring function.  c[i] is the current color of vertex i (value -1 means no color assigned).
    stack = []     # Our stack.  Item i is the set of vtxs we deleted color c[i] from (that is, the forward-nbrs which had c[i] in their lists).

    # At the top of this while loop, we have a partial coloring c.  We have either just colored vertex i-1 (in which case c has colored vertices 0 through i-1 and all of i's remaining colors are available to try), or we have backtracked from a partial coloring that didn't work (in which case we had colored i with something before, and now we need to try another of the remaining colors).
    while i < n:   # Go as long as there's something to be colored; failure of all colorings will throw an exception (empty stack) inside.

        # What colors have we not yet tried on i?  We've already tried all permissible colors of label at most c[i].
        remainingColors = {x for x in LA[i] if x>c[i]}  # (The comparison x>c[i] is why we use -1 for vertices with no assigned color.)

        backtrack = False

        # Check if any of the other uncolored vertices have no remaining colors.
        for v in range(i+1,n):
            if not LA[v]:
                backtrack = True

        # Case 1:  We can assign a new color to vertex i and we haven't discovered that farther vertices are empty.  So we do that.
        if remainingColors and not backtrack:
            c[i] = min(remainingColors)
            N = {x for x in G.neighbors(i) if x>i and x in CA[c[i]]}   # No need to bother with back-neighbors, since they have a color assigned.
            for v in N:
                LA[v].remove(c[i])    # Remove the new color from all forward-nbrs of i.
                CA[c[i]].remove(v)
            stack.append(N)
            i+=1

        # Case 2:  There are no colors left for some uncolored vertex, so we pop the stack and backtrack.
        else:
            try:
                oldN = stack.pop()
                for v in oldN:
                    LA[v].add(c[i-1])      # Add the attempted color back to the appropriate vertex lists.
                    CA[c[i-1]].add(v)
                c[i] = -1                # Mark i as uncolored.
                i-=1
            except IndexError:       # If the stack was empty, we've exhausted all possibilities.  Return an empty coloring.
                return []

    return c     # At this point, i == n.


n = 8
fin = open("geng_8.txt")

for line in fin:
    if line[0] == "<":
        continue
    graph6, lb, chromatic_num, ub = line.split()
    G = Graph(graph6)
    chromatic_num = int(chromatic_num)
    ub = int(ub)
    lb = int(lb)
    for i in xrange(lb, ub+1):
        f = [i]*n
        res = CNS(G,f)
        if res:
            at_num = i
            break
    if ub - at_num <= 1: # want gap of 2
        #print "Too small a gap", at_num, ub
        continue
    for i in xrange(chromatic_num, at_num+1):
        f = [i]*n
        res = fChoosable(G,f)
        if res[0]:
            print graph6, i, at_num, ub
            break

#family = ["E`~o", "GwC^~w", "I~?GW^~~o", "K~{?GKF@~~~}", "M~~w?CB?wF_^~~~~?", "O~~~w?@?WB_N?^?^~~~~}", "Q~~~~{??G@_F?N?N_Fw@~~~~~~o", "S~~~~~~???_B?F?F_Bw?~?F{?^~~~~~~w", "U~~~~~~~w??@?B?B_@w?^?B{?Nw?^w?^~~~~~~~o"]

#m = 2
#for string in family:
    #n = 2*m + 2
    #lb = int(m*(m+3)/(2*(m+1)))
    #G = Graph(string)
    #print "m =", m
    #print "n =", n
    #for i in xrange(lb,n):
        #f = [i]*n
        #res = CNS(G,f,verbosity = 8)
        #if res:
            #print i
            #print m+1
            #break
    #m += 1

#n = 7
#G = Graph("F?~v_")
#for i in xrange(2,n):
    #f = [i]*n
    #res = CNS(G,f,verbosity=8)
    #if res:
        #print i
        #break
