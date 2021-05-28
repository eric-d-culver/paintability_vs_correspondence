#!/usr/bin/python3

import itertools as it

# store a polynomial as a tuple and a dictionary
# tuple specifies the indices of the variables in this polynomial, and their order
# dictionary has keys of the tuples of exponents on those variables, and values of the coefficients
# example: x_4 + 5*x_4*x_5 + 17*x_5^2 would be stored as:
# (4,5), {(1,0): 1, (1,1): 5, (0,2): 17}

# returns true if x_i appears in poly
def var_appears(poly, i):
    if i in poly[0]:
        j = poly[0].index(i)
        for key in poly[1]:
            if key[j] > 0:
                return True
    else: 
        return False

# poly is a single polynomial
# returns single polynomial which is coefficient on x_i^e in the original
def coeff_extract_base(poly, i, e):
    if i in poly[0]:
        j = poly[0].index(i)
        new_vars = poly[0][:j] + poly[0][j+1:] # poly[0] without jth thing
        new_coeffs = {}
        for exp, coeff in poly[1]:
            if exp[j] == e:
                new_coeffs[exp[:j] + exp[j+1:]] = coeff
    else: # x_i is not present in poly
        if e == 0: # x_i^0 = 1, so want the whole thing
            return poly
        else:
            return 0

# Helper function
def equivalent_exponent(vars1, vars2, exp1):
    exp2 = [0]*len(vars2)
    for var, e in zip(vars1, exp1):
        if e != 0 and var not in vars2:
            return None # there is no equivalent
        elif e != 0:
            exp2[vars2.index(var)] = e
    return tuple(exp2)

# return expanded out product of polynomials
def expand_product(polys):
    pass

# return sum of polynomials
def sum_polys(polys):
    if len(polys) == 0:
        return ((), {(): 1}) # polynomial 1
    elif len(polys) == 1:
        return polys[0]
    else:
        first = polys[0]
        rest = sum_polys(polys[1:])
        res_coeffs = {}
        res_vars = tuple(set(first[0] + rest[0]))
        for exp1 in first[1]:
            exp2 = equivalent_exponent(first[0], rest[0], exp1)
            exp3 = equivalent_exponent(first[0], res_vars, exp1)
            if exp2 is not None and exp2 in rest[1]:
                res_coeffs[exp3] = first[1][exp1] + rest[1][exp2]
            else:
                res_coeffs[exp3] = first[1][exp1]
        for exp2 in rest[1]:
            exp3 = equivalent_exponent(rest[0], res_vars, exp2)
            if exp3 not in res_coeffs:
                res_coeffs[exp3] = rest[1][exp2]
        return (res_vars, res_coeffs)

# Helper function
# given e and n, yields all n-tuples of nonnegative integers which sum to e
def ntups_sum_to_e(e, n):
    for tup in it.product(range(e), repeat=n-1):
        last = e - sum(tup)
        if last >= 0:
            yield tup + (last,)

# polys is a list of polynomials which are being multiplied together
# this returns another list of polynomials being multiplied,
# which is the coefficient on x_i^e in the original polynomial product
# this partial expansion is the critical step of the algorithm, make it fast
def coeff_extract(polys, i, e):
    res = []
    relevant_factors = []
    for poly in polys:
        if var_appears(poly, i): # we need to consider this factor
            relevant_factors.append(poly)
        else: # if variable fails to appear, that factor carries over
            res.append(poly)
    # relevant_factors get expanded out to single factor, which is added to res
    for poss_exp in ntups_sum_to_e(e, len(relevant_factors)):
        pass # this is where the expanding out happens
    #res.append(new_factor)
    return res
