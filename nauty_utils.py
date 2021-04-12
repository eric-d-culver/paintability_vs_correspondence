#!/usr/bin/python3

# takes a filename, reads the file line by line interpreting each as a graph6 string, yields the adjacency list of each graph in turn
def graph6_file(filename):
    fin = open(filename, 'r')
    for line in fin:
        if line[0] == '<':
            continue
        yield graph6_to_dict(line[:-1])

# takes a graph6 format string representing a graph and outputs a dictionary which is the adjacency list of the graph (i.e., res[v] is a list of the neighbors of v)
def graph6_to_dict(g6_str):
    if len(g6_str) == 0:
        return # error
    if len(g6_str) == 1 and ord(g6_str[0]) == 63:
        return dict() # pointless graph
    if len(g6_str) == 1 and ord(g6_str[0]) == 64:
        return {0: []} # K1
    if len(g6_str) == 1:
        return # error 
    magic1 = ord(g6_str[0])
    magic2 = ord(g6_str[1])
    breakpt = 0
    if magic1 != 126:
        breakpt = 1
    elif magic2 != 126:
        breakpt = 4
    else: # magic1 == 126 and magic2 == 126
        breakpt = 8
    n = N_inv(g6_str[:breakpt])
    x = R_inv(g6_str[breakpt:])
    x = x[:int(n*(n-1)/2)]
    res = {i:[] for i in range(n)}
    k = 0 # k is index into x
    for i in range(n):
        for j in range(i):
            # i ranges from 0 to n-1
            # j ranges from 0 to i-1
            if x[k] == 1:
                # add edge from i to j
                res[i].append(j)
                res[j].append(i)
            k += 1
    return res


# inverse of N function from graph6 documentation
# returns integer between 0 and 2^36-1
def N_inv(string):
    if len(string) == 0:
        return # error
    if len(string) == 1:
        return ord(string[0]) - 63
    if len(string) == 4:
        return binary_num(R_inv(string[1:]))
    if len(string) == 8:
        return binary_num(R_inv(string[2:]))
    return # error

# inverse of R function from graph6 documentation
# returns list of bits
def R_inv(string):
    nums = [ord(c) - 63 for c in string]
    return sum(map(lambda a: pad_with_zeroes(binary_list(a), 6), nums), [])

# convert big-endian list of bits to number
def binary_num(lst):
    res = 0
    for i in lst:
        res += i
        res *= 2
    return res

# convert number to list of bits (big end first)
def binary_list(num):
    res = []
    while num != 0:
        res = [num % 2, *res]
        num = num // 2
    return res

# pad with zeroes at beginning
def pad_with_zeroes(bin_list, length):
    while len(bin_list) < length:
        bin_list = [0, *bin_list]
    return bin_list

#print(graph6_to_dict("DQ{"))

#for graph in graph6_file("geng_5.txt"):
    #print(graph)
