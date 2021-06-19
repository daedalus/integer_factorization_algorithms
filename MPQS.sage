#!/usr/bin/env python3
# Multiple polynomial Quadratic sieve
# Author Dario Clavijo 2021
# License: GPLv3

"""
The quadratic sieve algorithm (QS) is an integer factorization algorithm and, in practice, 
the second fastest method known (after the general number field sieve). 
It is still the fastest for integers under 100 decimal digits or so, and is considerably simpler than the number field sieve. 
It is a general-purpose factorization algorithm, meaning that its running time depends solely on the size of the integer to be factored,
and not on special structure or properties. 
It was invented by Carl Pomerance in 1981 as an improvement to Schroeppel's linear sieve
https://en.wikipedia.org/wiki/Quadratic_sieve
"""

import os
import time
import sys
from gmpy2 import gcd, gcdext, sqrt, isqrt, is_prime, next_prime, log2, log10, log, legendre, powmod, invert
from sage.parallel.multiprocessing_sage import parallel_iter
from multiprocessing import cpu_count, Pool, Manager
from itertools import repeat
import humanfriendly
from copy import *


def prebuilt_params(bits):
    """
    Bounds estimation
    borrowed from msieve/mpqs.c
    """
    if bits <= 64:
        return [100, 40, 1 * 65536]
    if bits <= 128: 
        return [450, 40, 1 * 65536]
    if bits <= 183:
        return [2000, 40, 1 * 65536]
    if bits <= 200: 
        return [3000, 50, 1 * 65536]
    if bits <= 212: 
        return [5400, 50, 3 * 65536]
    if bits <= 233:
        return [10000, 100, 3 * 65536]
    if bits <= 249:
        return [27000, 100, 3 * 65536]
    if bits <= 266:
        return [50000, 100, 3 * 65536]
    if bits <= 283:
        return [55000, 80, 3 * 65536]
    if bits <= 298:
        return [60000, 80, 9 * 65536]
    if bits <= 315:
        return [80000, 150, 9 * 65536]
    if bits <= 332:
        return [100000, 150, 9 * 65536]
    if bits <= 348:
        return [140000, 150, 9 * 65536]
    if bits <= 363:
        return [210000, 150, 13 * 65536]
    if bits <= 379:
        return [300000, 150, 17 * 65536]
    if bits <= 395:
        return [400000, 150, 21 * 65536]
    if bits <= 415:
        return [500000, 150, 25 * 65536] # beyond this point you're crazy 
    if bits <= 440:
        return [700000, 150, 33 * 65536]
    if bits <= 465:
        return [900000, 150, 50 * 65536]
    if bits <= 490:
        return [1100000, 150, 75 * 65536]
    if bits <= 512:
        return [1300000, 150, 100 * 65536]
    return [1300000, 150, 100 * 65536]


def choose_multiplier(n, prime_list):
    """
    Code borrowed from msieve/mpqs.c
    """
    mult_list = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38]
    mult_list += [39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73]

    MAX_MP_WORDS = 30
    NUM_MULTIPLIERS=len(mult_list)
    NUM_TEST_PRIMES = 300
    i =  j = 0
    best_score = 0.0
    best_mult = 0
    scores = [0.0] * NUM_MULTIPLIERS
    num_multipliers = 0
    M_LN2 = log(2)
    num_primes = len(prime_list)

    """ measure the contribution of 2 as a factor of sieve
       values. The multiplier itself must also be taken into
       account in the score. scores[i] is the correction that
       is implicitly applied to the size of sieve values for
       multiplier i; a negative score makes sieve values 
       smaller, and so is better """

    for i in range(0, NUM_MULTIPLIERS): 
        curr_mult = mult_list[i]
        knmod8 = (curr_mult * n) % 8
        logmult = log(curr_mult)
        # only consider multipliers k such than
        #   k*n will not overflow an mp_t */

        scores[i] = 0.5 * logmult;
        if knmod8 == 1:
            scores[i] -= 2 * M_LN2;
            break
        if knmod8 == 5:
            scores[i] -= M_LN2;
            break;
        if knmod8 == 3 or knmod8 == 7:
            scores[i] -= 0.5 * M_LN2;
            break;
        # even multipliers start with a handicap
    num_multipliers = i;

    # for the rest of the small factor base primes 

    for i in range (1, num_primes):
        prime = int(prime_list[i])
        contrib = log(prime) / (prime - 1);
        modp = n % prime

        for j in range(0, num_multipliers):
            curr_mult = mult_list[j];
            knmodp = (curr_mult * modp) % prime
            # if prime i is actually in the factor base
            #   for k * n ... */
            if (knmodp == 0 or legendre(knmodp, prime) == 1):

                """ ...add its contribution. A prime p con-
                   tributes log(p) to 1 in p sieve values, plus
                   log(p) to 1 in p^2 sieve values, etc. The
                   average contribution of all multiples of p 
                   to a random sieve value is thus
                   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
                   = (log(p) / p) * 1 / (1 - (1/p)) 
                   = log(p) / (p-1)
                   This contribution occurs once for each
                   square root used for sieving. There are two
                   roots for each factor base prime, unless
                   the prime divides k*n. In that case there 
                   is only one root """

                if (knmodp == 0):
                    scores[j] -= contrib
                else:
                    scores[j] -= 2 * contrib

    # use the multiplier that generates the best score 

    best_score = 1000.0
    best_mult = 1

    #print(scores)
   
    for i in range(0, num_multipliers):
        score = scores[i];
        if (score < best_score):
            best_score = score
            best_mult = mult_list[i]
    return best_mult;


def mod_sqrt(n, p):
    """ 
    Tonelli shanks algorithm 
    """
    a = n % p
    if p % 4 == 3:
        return pow(a, (p+1) >> 2, p)
    elif p % 8 == 5:
        v = pow(a << 1, (p-5) >> 3, p)
        i = ((a*v*v << 1) % p) - 1
        return (a*v*i) % p
    elif p % 8 == 1:  # Shank's method
        q, e = p-1, 0
        while q & 1 == 0:
            e += 1
            q >>= 1
        n = 2
        while legendre(n, p) != -1:
            n += 1
        w, x, y, r = pow(a, q, p), pow(a, (q+1) >> 1, p), pow(n, q, p), e
        while True:
            if w == 1:
                return x
            v, k = w, 0
            while v != 1 and k+1 < r:
                v = (v*v) % p
                k += 1
            if k == 0:
                return x
            d = pow(y, 1 << (r-k-1), p)
            x, y = (x*d) % p, (d*d) % p
            w, r = (w*y) % p, k
    else:
        return a  # p == 2


def trial_division(n, P):
    """ 
    A simple trial division, factors are given by P-list. 
    """
    a = []
    r = n
    l = len(P)
    i = 0
    for p in P:
        if r >= p:
            if r % p == 0:
                pw = 0
                while r % p == 0:
                    pw += 1
                    r //= p
                a.append((int(p),int(pw)))
        else:
            break
    return a,r,n


def merge_powers(ppws):
    """
    Merge powers such that: (a^b) * (a^c) == a^(b+c)
    """
    d = {}
    for p,pw in ppws:
        if p not in d:
            d[p] = pw
        else:
            d[p] += pw
    tmp2 = []
    for p in d:
       tmp2.append((p,d[p]))
    return tmp2


def filter_out_even_powers(ppws):
    """ 
    Same as merge but Filter out even powers. 
    """
    d = {}
    for p,pw in ppws:
        if p not in d:
            d[p] = pw
        else:
            d[p] += pw
    tmp2 = []
    for p in d:
        if d[p] & 1 != 0: # only keep odd powers
            tmp2.append(p)
    return tmp2


def trial_division_minus_even_powers(n, P):
    """
    Factor a composite n returning only odd power primes.
    """
    a, r, n = trial_division(n, P)
    a = filter_out_even_powers(a)
    return [a,r,n]


def is_smooth(x, P):
    """ 
    check if n is B-smooth. 
    """
    y = x
    for p in P:
        while y % p ==0:
            y //= p
    return abs(y) == 1 


def is_power(n, primes):
    """
    Given a interger n and a set of primes it checks if n is a perfect power.
    """
    ispow = False
    c = 0
    for p in primes:
        a = log(n) / log(p)
        b = int(a)
        if a == b:
            ispow = True
            c = b
    return ispow


def is_power_logprime(n, log_primes):
    """
    Given the precomputed logs of a set of primes 
    we can test if a number n is a perfect power of that prime.
    """
    ispow = False
    for p in log_primes:
        a = log(n) / p
        b = int(a)
        if a == b:
            iwpow = True
    return ispow


def minifactor4(x, P, smooth_base):
    """ 
    Minifactor. 
    """
    x = abs(x)
    a1 = []
    smooth = gcd(x, smooth_base)
    if smooth > 1:
        not_smooth = x // smooth
        a1, r1, n1 = trial_division(smooth, P)
        a2, r2, n2 = trial_division(not_smooth, P)
    else:
        a2, r2, n2 = trial_division(x, P) 
    return merge_powers(a1 + a2), r2, x


class Poly:
    """ 
    Quadratic polynomial helper class: 
    type Ax + 2Bx + C 
    """
    def __init__(self, n, P, x_max, search = 0, verbose = None, logp = None, A = None, B = None, C = None):
        self.n = int(n)
        self.P = P
        self.logp = logp
        self.x_max = x_max
        self.search = search
        self.verbose = verbose
        #self.counter = 0
        self.early_factors = []

        if A == None and B == None and C == None:
            self.create()
        else:
            print("given poly parameters: %d %d %d" % (A,B,C))
            self.A = A
            self.B = B
            self.C = C
            root_2n = isqrt(2 * self.n)
            root_A = next_prime(isqrt(root_2n // x_max))
            self.root_A = root_A

        self.solve_for_x()
 
 
    def create(self):
        """ 
        Construct our polynomial. 
        Code borrowed from 
        https://github.com/elliptic-shiho/primefac-fork/blob/master/_primefac/_factor_algo/_mpqs.py
        """
        n = self.n
        x_max = self.x_max
        root_2n = isqrt(2*n)
        root_A = next_prime(isqrt(root_2n // x_max))
        self.root_A = root_A
        s=0
        leg = 2
        while True:
            leg = legendre(n, root_A)
            if leg == 1:
                if s > self.search:
                    break
            elif leg == 0:
                self.early_factors.append(root_A)
                if s > self.search:
                    break
                #break
            s += 1
            root_A = next_prime(root_A + s)
        self.s = s
        A = int(pow(root_A, 2))
        # solve for an adequate B. B*B is a quadratic residue mod n,
        # such that B*B-A*C = n. this is unsolvable if n is not a
        # quadratic residue mod sqrt(A)
        b = int(mod_sqrt(n, root_A))
        B = int(b + (n - pow(b, 2)) * powmod((b * 2), (root_A-2) ,root_A)) % A
        C = int((pow(B, 2) - n) // A)        # B*B-A*C = n <=> C = (B*B-n)//A
        assert (B * B) - (A * C) == n
        assert C == ((B * B) - n) // A
        self.A = A
        self.B = B
        self.C = C
        if self.verbose:
            m = "New poly: f(x) = %dx^2 + %dx + %d\n" % (A,B,C)
            sys.stderr.write(m.replace("+ -","-"))
        return A, B, C


    def solve_for_x(self):
        """ 
        Get the minima of the polynomial where y=f(x) is a guaranteed relation. 
        Code borrowed https://github.com/cramppet/quadratic-sieve/blob/master/QS.py
        """        
        A = self.A
        B = self.B
        C = self.C
        n = self.n
        p = self.root_A
        if A > 1:
            g = 1
            while g == 1:
                p = next_prime(p)
                ainv = 1
                if A != 1:
                    g, inv, _ = gcdext(A, p)
                    if g != 1:
                        r1 = mod_sqrt(C, p)
                        r2 = (-1 * r1) % p
                        start1 = (ainv * (r1 - B)) % p
                        start2 = (ainv * (r2 - B)) % p
                        self.X = [int(start1), int(start2)]
                        self.minima = sum(self.X)//2
        elif A == 1:
            i = isqrt(abs(C))
            self.X = [-i,i]
            self.minima = 0


    def eval(self, x):
        """ 
        Eval the poly: y=f(x), return y and radical (Ax+B).
        """
        A = self.A
        B = self.B
        C = self.C
        Ax = A * x
        Rad = Ax + B 
        y = (Rad + B) * x + C
        return y, Rad, x # same as (Ax + 2B)x + C
      

    def __repr__(self):
        """
        Return the string representation of the constructed polynomial
        """
        m = "F(X) = %d X ^ 2 + %d X + %d with minima: %d" % (self.A, self.B, self.C, self.minima) 
        m = m.replace("+ -","- ")
        return m


    def __add__(self, other):
        """
        Adds one polynomial to other
        """
        if isinstance(other, Poly):
            if self.n == other.n and self.x_max == olther.x_max:
                A, B, C = self.A + other.A, self.B + other.B, self.C + other.C
                return Poly(self.n, self.x_max, search = None, A = A, B = B, C = C)
        else:
            return NotImplemented


    def __sub__(self, other):
        """
        Substracts one polynomial to other
        """
        if isinstance(other, Poly):
            if self.n == other.n and self.x_max == olther.x_max:
                A, B, C = self.A - other.A, self.B - other.B, self.C - other.C
                return Poly(self.n, self.x_max, search = None, A = A, B = B, C = C)
        else:
            return NotImplemented


    def __mul__(self, other):
        """
        Scalar multiply: scalar*(F(x))
        """
        if isinstance(other, int):
            A, B, C = self.A * other, self.B * other, self.C * other
            return Poly(self.n, self.x_max, search = None, A = A, B = B, C = C)
        else:
            return NotImplemented


    def __hash__(self):
        """
        Hashes unique values of the polynomial to construct an python internal id.
        """
        h = hash("%d-%d-%d" % (self.A,self.B,self.C))
        #print("hash: %s" % h)
        return h


    def __eq__(self, other):
        """
        Internal python facility needed to check if two polynomials are the same.
        """
        if isinstance(other, Poly):
            return ((self.A == other.A) and (self.B == other.B) and (self.C == other.C))
        else:
            return NotImplemented


def compute_logs_y_sums(start, stop, P, log_primes, last_logp_table = []):
     """
     Compute log table
     """
     logs_y = copy(last_logp_table)
     logs_y += [0] * (stop-start)
     for j in range(len(P)):
         p = P[j]
         log_p = log_primes[j]
         for i in range(0 , stop , p):
             if i >= start:
                 logs_y[i] += log_p
     return logs_y


def relations_find(taskid, N, start, stop, P, min_log_primes, log_primes, logs_y, smooth_base, Rels, merged_count, required_relations, cycleFactors, thresh, tasks, polycounts, polys = None):
    """ 
    Relations search funcion 
    """
    #pid = os.getpid()
    pid = taskid
    tasks.value += 1
    sys.stderr.write("[%d] relations_find: range(%d, %d), interval: %d sieving start, tresh: %f, pb: %d.\n" % (pid,start, stop, (stop-start), thresh,len(polys)))
    pre_log_filter = True
    merge = True
    proc = noproc = 0
    Found_Rels = []
    st = time.time()
    ltd = st
    m = 1000
    msg = ""
    partials = {}
    rels_found = 0

    for poly in polys:
        m = poly.minima
        Diffs = [poly.eval(abs(x)) for x in range(start + m, stop +m)]
        ld = len(Diffs)
        A = poly.A
        for i in range(ld):
            if len(Rels) > required_relations or len(cycleFactors) > 0:
                polycounts[poly] = (start,stop,rels_found)
                break
            y, Rad, x = Diffs[i]
            y = abs(y)

            if pre_log_filter:
                candidate = (0 < logs_y[i] <= thresh and not is_prime(y))
            else:
                candidate = (y > P[0] and not is_prime(y) and not is_square(y) and not is_power_logprime(y, min_log_primes))

            if candidate:
                proc += 1
                f = minifactor4(y, P, smooth_base)
                if f != None:            
                    if f[1] == 1:  # Found a relation 
                        rels_found += 1
                        filtered = filter_out_even_powers(f[0])
                        rel = [filtered, y, Rad, A]
                        if rel not in Rels:
                            Rels.append(rel)
                    elif merge and (f[1] in partials): # Found a partial, try to find a cycle
                        a = partials[f[1]]
                        p = filter_out_even_powers(f[0] + a[0])
                        Ahs = A * a[3]
                        lhs = Rad * a[2]
                        rhs = y * a[1]
                        LHS = Ahs + lhs
                        g = gcd(isqrt(LHS) - rhs, N)
                        if N > g > 1:
                            cycleFactors.append([g , N // g])
                            sys.stderr.write("Found cycle with partial\n")
                        else: # Cycle not found, merge
                            if merge:
                                rel = [p, rhs, lhs, Ahs]
                                if rel not in Rels:
                                    Rels.append(rel)
                                    rels_found += 1
                                    #with merged_count.value.get_lock():
                                    merged_count.value += 1
                                    del partials[f[1]]
                    else:
                        if merge: # Didnt merge, store the partial for later merging
                            partials[f[1]] = [f[0], y, Rad, A]
            else:
                noproc += 1

            if i % m == 0:
                lRels = len(Rels)
                lt = time.time()
                td = lt - ltd
                ltd = lt
                if i > 0:
                    if lRels > 0:
                        eta = int(td * sqrt((ld / m)**2 + (required_relations / lRels)**2))
                    else:
                        eta = int(td * (ld / m))
                    etas = humanfriendly.format_timespan(eta)
                else:
                    etas = "No ETA estimated"
                tds = humanfriendly.format_timespan(td)
                msg = "[%d] relations_find: range (%d, %d , %d), found: %d of %d, merged: %d, proc: %d, noproc: %d, iter_elapsed: %s, eta: %s, alive tasks: %d.\n" 
                msg = msg % (pid, start, i, ld, lRels, required_relations, merged_count.value, proc, noproc , tds, etas, tasks.value)
                sys.stderr.write(msg)
    
        polycounts[poly] = (start,stop,rels_found)

    td = time.time() - st
    tds = humanfriendly.format_timespan(td)
    msg = "[%d] relations_find: Ended range(%d, %d, %d), found: %d of %d, merged: %d, proc: %d, noproc: %d, time elapsed: %s, alive tasks: %d.\n" 
    msg = msg % (pid, start,i,ld,len(Rels),required_relations, merged_count.value, proc, noproc, tds, tasks.value)
    sys.stderr.write(msg)
    tasks.value -= 1

    return Found_Rels
    

def transpose(A):
    """
    Transpose matrix so columns become rows
    """
    new_A = []
    for i in range(len(A[0])):
        new_row = []
        for row in A:
            new_row.append(row[i])
        new_A.append(new_row)
    return(new_A)


def Gaussian_elimination_GF2(A):
  """
  Gaussian elimination over GF2
  https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf, page2
  """
  h = len(A)
  m = len(A[0])
  marks = [False] * h
  for j in range(0,m):
    for i in range(0,h):
      if A[i][j] == 1:
          marks[i] = True
          for k in range(j+1,m):
            if A[i][k] == 1:
              A[i][k] = (A[i][j] + A[i][k]) % 2 
          break
  return marks, A


def create_matrix(Rels, P):
    """
    Create a matrix with the relations.
    """
    M1 = []
    for i in range(0,len(Rels)):
        row = []
        for j in range(0,len(P)): 
            if P[j] in Rels[i][0]:
                row.append(1)
            else:
                row.append(0)
        M1.append(row)
    return M1


def left_nullspace(A):
    """
    Compute left null space:
    x such that xTA = 0T
    """
    A = transpose(A)
    marks, A = Gaussian_elimination_GF2(A)
    #marks, A = gauss(A)

    B = []
    for row_index in range(0,len(A)):
        if marks[row_index] == False:
            B.append(A[row_index])
    return B


def linear_algebra(Rels, P):
    """ 
    Linear algebra, it generates a matrix in GF(2) then computes it's left null-space. 
    """
    M = matrix(GF(2), len(Rels), len(P), lambda i, j:P[j] in Rels[i][0])
    return M.left_kernel().basis()


def process_basis_vectors(N, basis, Rels, multiplier = 1):
    """ 
    Process each basis vector, construct (a^2)-(b^2) mod n congruence. 
    """
    for K in basis:
        lhs = rhs = Ahs = 1
        I = [f for f, k in zip(Rels, K) if k == 1]
        if len(I) > 0:
            for i in I:
                lhs *= i[1] # left-hand side
                rhs *= i[2] # rigth-hand side
                Ahs *= i[3] # A-term in poly
            LHS = Ahs * lhs
            g = gcd(isqrt(abs(LHS))-rhs,N)
            if N > g > 1:     
                factors = [int(g), int(N//g)]
                if multiplier > 1:
                    tmp2 = []
                    for factor in factors:
                        ignoreme  = gcd(factor,multiplier)
                        if r > ignoreme > 1:
                            tmp2.append(factor//ignoreme)
                    return tmp2
                else:
                    return factors


def find_primebase(n, bound):
    """ 
    Finds the base prime for given n and bound. 
    https://github.com/elliptic-shiho/primefac-fork/blob/master/_primefac/_factor_algo/_mpqs.py
    Same as the sieve of erathostenes.
    """
    print(n,bound)
    primes, mod_root, log_p, num_prime, p = [], {}, {}, 0, 3
    while num_prime <= bound:
        leg = legendre(n % p, p)
        if leg == 1:
            primes += [p]
            #mod_root[p] = mod_sqrt(n ,p)
            log_p[p] = log10(p)
            num_prime += 1
        elif leg == 0:
            return [p], None
        p = next_prime(p)
    return primes, log_p


def recalculate_min_prime_thresh(thresh, Prime_base, log_p):
    """ 
    Recalculates the min efective prime base bound. 
    code borrowed from:
    https://github.com/elliptic-shiho/primefac-fork/blob/master/_primefac/_factor_algo/_mpqs.py
    """
    min_prime = int(thresh * 3)
    fudge = sum(log_p[i] for i, p in enumerate(Prime_base) if p < min_prime)
    fudge = fudge // 4
    sys.stderr.write("min_prime: %d, fudge: %f, thresh: %f\n" % (min_prime,fudge, thresh))
    thresh -= fudge
    return min_prime, thresh, fudge


def generate_polys(N, Prime_base, x_max, needed, min_search = 0, polys=[]):
    """ 
    It searchs for distinct needed polys congruent to n. 
    """
    n = min_search
    cpolys = 0
    early_factors = []
    while cpolys <= needed:
        pol = Poly(N, Prime_base, x_max, search = n, verbose=False)
        if len(pol.early_factors) > 0:
            for early_factor in pol.early_factors:
                if early_factor not in early_factors:
                    early_factors.append(early_factor)
        else:
            if pol not in polys and pol.minima > 0:
                m = repr(pol) 
                sys.stderr.write("New Poly %d: %s\n" % (cpolys,m))
                polys.append(pol)
                cpolys += 1
        n += 1
    return polys, early_factors


def getbestpolys(polys, polycounts, T):
    """
    Experimental function, break things.
    """
    tmp_polys = copy(polys)
    new_polys = []
    print("have polys: %d" % len(tmp_polys))
    print(polys)
    j = -1
    while len(new_polys) < T:
        max_poly_count = 0
        for i in range(0,len(tmp_polys)):
            poly = tmp_polys[i]
            if poly in polycounts:
                poly_count = polycounts[poly][2]
                if poly_count > max_poly_count:
                    max_poly_count = poly_count
                    best_poly = poly                
                    j = i

        #if -1 < j < len(tmp_polys):
        if -1 < j < len(tmp_polys):
            print("new best poly: %s, j: %d, len: %d, count: %d" % (repr(poly),j,len(tmp_polys),max_poly_count))
            new_polys.append(best_poly)
            del tmp_polys[j]
        else:
            print(tmp_polys)
            break
    sys.stderr.write("Requested best_polys: %d, got: %d\n" % (T,len(new_polys)))
    return new_polys


def poly_stats(polys, polycounts):
    """
    Prints polynomials relations stats
    """
    for poly in polys:
        if poly in polycounts:
            poly_count = polycounts[poly]
            sys.stderr.write("For polynomial: %s the relations count is: %s\n" % (repr(poly),poly_count))


def unique(List, newlist):
    """
    Creates a unique list of elements
    """
    for element in List:
        if element not in tmp:
            newlist.append(element)
        else:
            print("repeated")
    return newlist
    

def _MPQS(N, verbose=True, M = 5):
    """ 
    Main MPQS function. 
    """
    if is_prime(N):
        return [N], 0 

    bN, lN = int(log2(N)), len(str(N))
    i2N = isqrt(N)

    if pow(i2N,2) == N:
        return [i2N, i2N],0

    i2Np1 = i2N + 1 
    root_2n = isqrt(2*N)
    Rels = []

    T = cpu_count() * M
    #T *= M

    B2, _ , B1 = prebuilt_params(log2(N))
    
    x_max = B1
    m_val = (x_max * root_2n) >> 1
    thresh = log10(m_val) * 0.735
    min_prime = int(thresh * 3)

    Prime_base, log_p = find_primebase(N, B2 + min_prime)

    multiplier = choose_multiplier(N, Prime_base)
    Nm = multiplier * N
    
    required_relations = len(Prime_base)

    if required_relations == 1:
        sys.stderr.write("Found small factor: %d\n" % Prime_base[0])
        r, polycount = _MPQS(N // Prime_base[0])
        return Prime_base + r, polycount

    sys.stderr.write("Factoring N: %d, bits: %d, digits: %d, B1: %d, B2: %d\n" % (N,bN,lN,B1,B2))
    sys.stderr.write("Multiplier is: %d\n" % multiplier)

    start = 0
    stop = B1 # range to sieve
    
    manager = Manager()
    Rels = manager.list() # placeholder list for relations shareable between child processes.
    newlist = manager.list()
    cycleFactors = manager.list()
    polycounts = manager.dict()

    merged_count = manager.Value("i", 0)
    tasks = manager.Value("i", 0)

    min_poly_search = 0
    polys = []

    last_lRels = 0
    force_new_polys = True
    use_best = True

    required_relations_ratio = 1.05
    required_relations = int(required_relations * required_relations_ratio)

    logs_y = []

    while True:
        # trim primes, recalc min

        sys.stderr.write("Generating primebase and logprime table...\n")
        min_log_primes = [log(p) for p in Prime_base if p <= min_prime]
        filtered_Prime_base = [p for p in Prime_base if p > min_prime]
        log_primes = [log(p) for p in Prime_base]
        filtered_log_primes = [log(p) for p in filtered_Prime_base]
      
        if stop > len(logs_y):
            logs_y = compute_logs_y_sums(start, stop, filtered_Prime_base, filtered_log_primes, logs_y)
            required_relations = int(required_relations * required_relations_ratio)

        sys.stderr.write("Need %d relations\n" % (required_relations))
        
        smooth_base = prod(Prime_base)
        min_prime, thresh, fudge = recalculate_min_prime_thresh(thresh, filtered_Prime_base, log_p)

        t1 = time.time()
        sys.stderr.write("Data collection with %d threads...\n" % T)

        lRels = len(Rels)

        need_more_polys = not (lRels > last_lRels) or force_new_polys
        need_more_polys = False

        last_lRels = lRels

        inputs = [] 
        taskid = 1

        # generate tasks parameters
        if not need_more_polys:
            if len(polys) > 0 and use_best:
                sys.stderr.write("Using the best polys to get more relations...\n")
                polys_new = getbestpolys(polys, polycounts, T)
                polys_new += polys[0 - T:len(polys)]
            else:
                #polys_new = generate_polys(Nm, Prime_base, x_max, T, min_poly_search, polys)
                sys.stderr.write("Need more polys...\n")
                polys_new, early_factors = generate_polys(Nm, Prime_base, x_max, T, min_poly_search, polys) # generate n distinct polys one for each cpu core.
                if len(early_factors) > 0:
                    tmp = 1
                    small = []
                    for early_factor in early_factors:
                        g = gcd(early_factor, N)
                        if N > g > 1:
                            tmp *= g
                            small.append(g)
                    if tmp > 1:
                        return small + MPQS(N // tmp), 0
                min_poly_search += T
        else:
            polys_new = polys[0 - (T * M):len(polys)]
            
        pb = M
        for i in range(0, len(polys_new), pb):
            polys_batch = polys_new[i : i + pb]

            inputs += [(taskid, Nm, start, stop, filtered_Prime_base, min_log_primes, log_primes, logs_y, smooth_base, Rels, merged_count, required_relations, cycleFactors, thresh, tasks, polycounts, polys_batch)]
            taskid += 1
            for poly in polys_batch:
                if poly not in polys:
                    poly.append(poly) 

        # deploy tasks to every cpu core.
        pols = []
        workpool = Pool(T)
        with workpool:
            R = workpool.starmap(relations_find, inputs)  
        workpool.close()
        workpool.join() 
        
        t2 = time.time()
        sys.stderr.write("Done in: %f secs.\n" % (t2-t1))
        #sys.stderr.write("Found %d rels with %d base primes.\nSorting..." % (len(Rels),len(Prime_base)))
        sys.stderr.write("Found %d rels with %d base primes.\n" % (len(Rels),len(Prime_base)))


        if len(cycleFactors) > 0:
            return cycleFactors, len(polys)

        t3 = time.time()

        #sys.stderr.write("Done in: %f secs.\n" % (t3-t2))
        #sys.stderr.write("Unique rels after sort: %d \n" % (len(Rels)))
        
        # when needed relations is reached proceed to do linear algebra
        if len(Rels) > required_relations:
            sys.stderr.write("Found %d enough relations of %d needed relations...\n" % (len(Rels),required_relations))
            sys.stderr.write("Matrix creation...")
            
            basis = linear_algebra(Rels, Prime_base) # calculate left nullspace of a GF(2) matrix.

            t4 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t4-t3))
            sys.stderr.write("Matrix reduction...")
            result = process_basis_vectors(Nm, basis, Rels, multiplier)
            t5 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t5-t4))
             
            poly_stats(polys, polycounts)

            if result != None:
                return result, len(polys)

        if not need_more_polys:
            start = stop
            stop += B1
        else:
            min_poly_search += T

        sys.stderr.write("Need to sieve %d differences more...\nNew sieving range %d:%d\n" % (required_relations - len(Rels),start,stop))
      

def MPQS(N):
    """ 
    Iterative version of MPQS. 
    """
    result, polys = _MPQS(N)
    polycount = polys

    if result != None:
        R= []
        for r in result:
            if is_prime(r):
                sys.stderr.write("Found factor: %d\n" % r)
                R.append(r)
            else:
                r, polys = _MPQS(r)
                R+=r
                polycount += polys
        return R, polycount


if __name__ == "__main__":
    N = Integer(sys.argv[1])
    ts = time.time()
    r, polycount = MPQS(N)
    Fp = trial_division(N,sorted(r))[0] 
    tmp = []
    for f,p in Fp:
        tmp += ["%d ^ %d" % (f,p)]
    msg = "The factorization for N = %d is: %s" % (N," * ".join(tmp))
    print(msg.replace(" ^ 1",""))
    print("With %d distinct generator polys" % polycount)
    td = time.time() - ts
    tds = humanfriendly.format_timespan(td) 
    sys.stderr.write("All done in: %s.\n" % (tds))
