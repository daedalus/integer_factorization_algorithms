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


def is_power_logprime(n, min_log_primes):
    ispow = False
    for p in min_log_primes:
        a = log(n) / p
        b = int(a)
        if a == b:
            iwpow = True
    return ispow


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
    Filter out even powers. 
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
    def __init__(self, n, P, x_max, search = 0, verbose = None, logp = None):
        self.n = int(n)
        self.P = P
        self.logp = logp
        self.x_max = x_max
        self.search = search
        self.verbose = verbose
        #self.counter = 0
        self.early_factors = []
        self.create()
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
            root_A = next_prime(root_A + s)
            #root_A = next_prime(root_A)

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


    def eval(self, x):
        """ 
        Eval the poly: y=f(x), return y and radical (Ax+B).
        """
        A = self.A
        B = self.B
        C = self.C
        Ax = A * x
        Rad = Ax + B 
        return (Rad + B) * x + C, Rad # same as (Ax + 2B)x + C
      
    def __repr__(self):
        """
        Return the string representation of the constructed polynomial
        """
        m = "F(X) = %d X ^ 2 + %d X + %d with minima: %d" % (self.A, self.B, self.C, self.minima) 
        m = m.replace("+ -","- ")
        return m

    def __hash__(self):
        h = hash("%d-%d-%d" % (self.A,self.B,self.C))
        #print("hash: %s" % h)
        return h

    def __eq__(self, other):
        if isinstance(other, Poly):
            return ((self.A == other.A) and (self.B == other.B) and (self.C == other.C))
        else:
            return NotImplemented

def is_power(n, minprimes):
    ispow = False
    c = 0
    for p in minprimes:
        a = log(n) / log(p)
        b = int(a)
        if a == b:
            ispow = True
            c = b
    return ispow


def is_power_logprime(n, min_log_primes):
    """
    Given the precomputed logs of a set of primes 
    we can test if a number n is a perfect power of that prime.
    """
    ispow = False
    for p in min_log_primes:
        a = log(n) / p
        b = int(a)
        if a == b:
            iwpow = True
    return ispow


def relations_find(taskid, N, start, stop, P, min_log_primes, smooth_base, Rels, merged_count, required_relations, cycleFactors, thresh, tasks, polycounts, poly = None):
    """ 
    Relations search funcion 
    """
    #pid = os.getpid()
    pid = taskid
    tasks.value += 1
    sys.stderr.write("[%d] relations_find: range(%d, %d), interval: %d sieving start\n" % (pid,start, stop, (stop-start)))
    m = poly.minima
    Diffs = [(poly.eval(abs(x)),abs(x)) for x in range(start + m, stop +m)]
    proc = noproc = 0
    Found_Rels = []
    A = poly.A
    st = time.time()
    ltd = st
    ld = len(Diffs)
    m = 1000
    msg = ""
    partials = {}
    i= 0
    rels_found = 0
    for i in range(ld):
        if len(Rels) > required_relations or len(cycleFactors) > 0:
            break
        yRad, x = Diffs[i]
        y, Rad = yRad
        y = abs(y)
        if y > P[0] and not is_prime(y) and not is_square(y) and not is_power_logprime(y, min_log_primes):
            proc += 1 
            f = minifactor4(y, P, smooth_base)
            if f != None:            
                if f[1] == 1:  # found a relation 
                    rels_found += 1
                    filtered = filter_out_even_powers(f[0])
                    rel = [filtered, y, Rad, A] 
                    #print(rel)
                    Rels.append(rel)
                elif f[1] in partials: # found a partial and try to merge
                    a = partials[f[1]]
                    p = filter_out_even_powers(f[0] + a[0])
                    Ahs = A*a[3]
                    lhs = Rad * a[2]
                    rhs = y * a[1]
                    LHS = Ahs + lhs
                    g = gcd(isqrt(LHS) - rhs, N)
                    if N > g > 1:
                        cycleFactors.append([g , N // g])
                        sys.stderr.write("Found cycle with partial\n")
                    else: # didnt merge so store the partial for later merging
                        Rels.append([p, rhs, lhs, Ahs])
                        rels_found += 1
                    #with merged_count.value.get_lock():
                    merged_count.value += 1
                else:
                    partials[f[1]] = [f[0], y, Rad, A]
        else:
            noproc += 1

        if i % m == 0:
            #print(repr(poly))
            #poly.counter = rels_found
            polycounts[poly] = rels_found
            lRels = len(Rels)
            lt = time.time()
            td = lt - ltd
            ltd = lt
            if i > 0:
                if lRels > 0:
                    #eta = td * (((ld / m) + (required_relations / lRels)) / 2)
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
    
    polycounts[poly] = rels_found
    td = time.time() - st
    tds = humanfriendly.format_timespan(td)
    msg = "[%d] relations_find: Ended range(%d, %d, %d), found: %d of %d, merged: %d, proc: %d, noproc: %d, time elapsed: %s, alive tasks: %d.\n" 
    msg = msg % (pid, start,i,ld,len(Rels),required_relations, merged_count.value, proc, noproc, tds, tasks.value)
    sys.stderr.write(msg)
    tasks.value -= 1
    #Rels += Found_Rels
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
              A[i][k] = (A[i][j] ^ A[i][k]) 
          break
  return marks, A

def left_nullspace(A):
    """
    Compute left null space:
    x such that xTA = 0T
    """
    A = transpose(A)
    marks, A = Gaussian_elimination_GF2(A)
    #marks, A = gauss(A)

    B = []
    for row in range(len(A)):
        if marks[row] == False:
            B.append(A[row])
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


#def generate_polys(N, Prime_base, x_max, needed, min_search = 0, polys_ABC=[]):
def generate_polys(N, Prime_base, x_max, needed, min_search = 0, polys=[]):

    """ 
    It searchs for distinct needed polys congruent to n. 
    """
    n = min_search
    cpolys = 0
    #polys = []
    early_factors = []
    while cpolys <= needed:
        pol = Poly(N, Prime_base, x_max, search = n, verbose=False)
        if len(pol.early_factors) > 0:
            for early_factor in pol.early_factors:
                if early_factor not in early_factors:
                    early_factors.append(early_factor)
        else:
            #pol_ABC = (pol.A,pol.B,pol.C)
            #if pol_ABC not in polys_ABC:
            if pol not in polys:
                #m = "Found poly: f(x) = %d X ^ 2 + %d X + %d, minima: %d\n" % (pol_ABC[0], pol_ABC[1], pol_ABC[2], pol.minima)
                m = repr(pol) 
                sys.stderr.write("New Poly %d: %s\n" % (cpolys,m))
                #polys_ABC.append(pol_ABC)
                polys.append(pol)
                cpolys += 1
        n += 1
    #return polys, polys_ABC, early_factors
    return polys, early_factors



def poly_stats(polys, polycounts):
    for poly in polys:
        if poly in polycounts:
            poly_count = polycounts[poly]
            sys.stderr.write("For poly: %s the count is: %d\n" % (repr(poly),poly_count))


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


def _MPQS(N, verbose=True, M = 1):
    """ 
    Main MPQS function. 
    """
    if is_prime(N):
        return [N]

    bN, lN = int(log2(N)), len(str(N))
    i2N = isqrt(N)

    if pow(i2N,2) == N:
        return [i2N, i2N]

    i2Np1 = i2N + 1 
    root_2n = isqrt(2*N)
    Rels = []

    T = cpu_count()

    B2, _ , B1 = prebuilt_params(log2(N))
    Prime_base, log_p = find_primebase(N, B2)
    #B1 //= M
    B2 = len(Prime_base)
  
    if B2 == 1:
        sys.stderr.write("Found small factor: %d\n" % Prime_base[0])
        return Prime_base + _MPQS(N // Prime_base[0]), 0

    multiplier = choose_multiplier(N, Prime_base)
    Nm = multiplier * N
    
    #x_max = B2 *60  # size of the sieve
    x_max = B1
    m_val = (x_max * root_2n) >> 1
    thresh = log10(m_val) * 0.735
    min_prime = int(thresh * 3)
    


    sys.stderr.write("Factoring N: %d, bits: %d, digits: %d, B1: %d, B2: %d\n" % (N,bN,lN,B1,B2))
    sys.stderr.write("Multiplier is: %d\n" % multiplier)

    start = 0
    stop = B1 # range to sieve
    
    manager = Manager()
    Rels = manager.list() # placeholder list for relations shareable between child processes.
    cycleFactors = manager.list()
    polycounts = manager.dict()

    merged_count = manager.Value("i", 0)
    tasks = manager.Value("i", 0)

    min_poly_search = 0
    #polys_ABC = []
    polys = []

    last_lRels = 0
    force_new_polys = False

    required_relations_ratio = 1.05
    required_relations = len(Prime_base)

    while True:
        # trim primes, recalc min
        required_relations = int(required_relations * required_relations_ratio)
        sys.stderr.write("Need %d relations\n" % (required_relations))
        min_log_primes = [log(p) for p in Prime_base if p <= min_prime]
        Prime_base = [p for p in Prime_base if p > min_prime]
        log_primes = [log(p) for p in Prime_base]
        smooth_base = prod(Prime_base)
        min_prime, thresh, fudge = recalculate_min_prime_thresh(thresh, Prime_base, log_p)

        t1 = time.time()
        sys.stderr.write("Data collection with %d threads...\n" % T)

        lRels = len(Rels)
        need_more_polys = not (lRels > last_lRels) or force_new_polys
        last_lRels = lRels

        if need_more_polys == True: 
            sys.stderr.write("Need more polys...\n") 
            #polys, polys_ABC, early_factors = generate_polys(Nm, Prime_base, x_max, T * 2, min_poly_search, polys_ABC) # generate n distinct polys one for each cpu core.
            polys, early_factors = generate_polys(Nm, Prime_base, x_max, T * 2, min_poly_search, polys) # generate n distinct polys one for each cpu core.

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
             
        inputs = [] 
        taskid = 1
        # generate tasks parameters
        for poly in polys[0 - (T * 2):len(polys)]:
            #print(poly.sums)
            #calc_sums
            inputs += [(taskid, Nm, start, stop, Prime_base, min_log_primes, smooth_base, Rels, merged_count, required_relations, cycleFactors, thresh, tasks, polycounts, poly)]
            taskid += 1
        #sys.exit(0)

        # deploy tasks to every cpu core.
        pols = []
        workpool = Pool(T)
        #print(T,len(inputs))
        #sys.exit(0)
        with workpool:
            R = workpool.starmap(relations_find, inputs)  
        workpool.close()
        workpool.join() 
        
        t2 = time.time()
        sys.stderr.write("Done in: %f secs.\n" % (t2-t1))
        sys.stderr.write("Found %d rels with %d base primes.\n" % (len(Rels),len(Prime_base)))

        if len(cycleFactors) > 0:
            return cycleFactors

        # when needed relations is reached proceed to do linear algebra
        if len(Rels) > required_relations:
            sys.stderr.write("Found %d enough relations of %d needed relations...\n" % (len(Rels),required_relations))
            sys.stderr.write("Matrix creation...")
            
            basis = linear_algebra(Rels, Prime_base) # calculate left nullspace of a GF(2) matrix.

            t3 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t3-t2))
            sys.stderr.write("Matrix reduction...\n")
            result = process_basis_vectors(Nm, basis, Rels, multiplier)
            t4 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t4-t3))
             
            poly_stats(polys, polycounts)

            if result != None:
                return result, len(polys)

        if not need_more_polys:
            start = stop
            stop += B1
        else:
            min_poly_search += (T * 2)

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
