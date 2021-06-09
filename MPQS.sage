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

import time
import sys
from gmpy2 import gcd, gcdext, isqrt, is_prime, next_prime, log2, log10, legendre, powmod, invert
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


def trial_division(n, P):
    """ 
    A simple trial division, factors are given by P-list. 
    """
    a = []
    r = n
    l = len(P)
    i = 0
    pi2 = pow(P[i],2)
    while (i < l) and (pi2 <= n):
        if r % P[i] == 0:
            pw = 0
            while r % P[i] == 0:
                pw += 1
                r //= P[i]
            a.append((int(P[i]),int(pw)))
        else:
            i += 1
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


def minifactor2(x, P, D = 1):
    """ 
    Tries to factor out common primes with gcd and filter out even power primes. 
    """
    tmp = []
    R = 1
    g = gcd(x, D) 
    if x > g > 1: # if there is a common factor try to factor primes from that factor
        gtmp = [g, x//g]
        for g1 in gtmp:
            if not is_prime(g1):
                i = isqrt(abs(g1))
                if i**2 != abs(g1): # g1 should not be square
                    a, r, n = trial_division(g1,P)
                    tmp += a
                    R *= r
                else:
                    if is_prime(i):
                        if i > P[-1] or i < P[0]:  
                            R *= i
                        else:
                        #a, r, n = [i], 1, g1
                            tmp += [(i,2)]
                            #R //= g1
                    else:
                        #R *= r
                        gtmp.append(i)    
            else:
                if g1 > P[-1] or g1 < P[0]:
                    R *= g1
                else:
                    tmp += [(g1,1)]
                    #R //= g1
        D//=g # reduce D
    else: # if not proceed as always
        a, r, n = trial_division(x,P)
        R *= r
        tmp += a
        D//=x # reduce D
    #print(tmp)
    if R == 1: # if remainder is 1 then we found a B-smooth number then relation.
        tmp = filter_out_even_powers(a)  
        return (tmp, R, x), D                 


def minifactor3(x, P, smooth_base):
    """ 
    Minifactor. 
    """
    x = abs(x)
    smooth = gcd(x, smooth_base)
    if smooth > 1:
        not_smooth = x // smooth
        a1, r1, n1 = trial_division(smooth, P)
        a2, r2, n2 = trial_division(not_smooth, P)
        if r2 == 1:
            a = a1 + a2
            a3 = filter_out_even_powers(a)
            return a3, r2, x


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

        
                         
def minifactor(x, P):
    """ 
    Minifactor algo, finds odd-power primes in composites. 
    """
    p = trial_division_minus_even_powers(x,P)
    if p[1] == 1: # if 1 x is B-smooth
        return p


class Poly:
    """ 
    Quadratic polynomial helper class: 
    type Ax + 2Bx + C 
    """
    def __init__(self, n, P, x_max, search = 0, verbose = None):
        self.n = int(n)
        self.P = P
        self.x_max = x_max
        self.search = search
        self.verbose = verbose
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
            root_A = next_prime(root_A)
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
        Get the mininmum x value of the polynomial where y=f(x) is a guaranteed relation. 
        Code borrowed https://github.com/cramppet/quadratic-sieve/blob/master/QS.py
        """
        
        A = self.A
        B = self.B
        C = self.C
        n = self.n
        start_vals = []
        #for p in self.P:
        p = self.root_A
        g = 1
        while g == 1:
        #for p in self.P:
            p = next_prime(p)
            ainv = 1
            if A != 1:
                 g, inv, _ = gcdext(A, p)
            if g != 1:
                r1 = mod_sqrt(C, p)
                r2 = (-1 * r1) % p
                start1 = (ainv * (r1 - B)) % p
                start2 = (ainv * (r2 - B)) % p
                start_vals.append([int(start1), int(start2)])
                self.start_vals = start_vals
        return start_vals


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
        

def relations_find(N, start, stop, P, smooth_base, Rels, merged_count, required_relations, pol = None):
    """ 
    Relations search funcion 
    """
    sys.stderr.write("relations_find: range(%d, %d), interval: %d sieving start\n" % (start, stop, (stop-start)))

    #if (stop-start) < 0:
    #    return [] 

    if (stop-start) > 0:
        I = [abs(x) for x in range(start, stop) if not is_prime(x) and not is_square(x)]
    else: # this may not happen
        I = [abs(x) for x in range(start, stop, -1) if not is_prime(x) and not is_square(x)]

    #D = reduce(lambda x, y: x * y, I)
    Diffs = [(pol.eval(x),x) for x in I]
    Found_Rels = []
    A = pol.A
    #B = pol.B
    #C = pol.C
    st = time.time
    ltd = st
    ld = len(Diffs)
    m = 1000
    msg = ""
    partials = {}

    for i in range(ld):
        if len(Rels) > required_relations:
            break
        yRad, x = Diffs[i]
        y, Rad = yRad
        f = minifactor4(y, P, smooth_base)
        if 1:
            if f != None:
                filtered = (filter_out_even_powers(f[0]),f[1],y)
                if f[1] == 1:    
                    Rels.append([filtered, y, Rad, A])
                elif f[1] in partials:
                    a = partials[f[1]]
                    p = filter_out_even_powers(f[0] + a[0])
                    Rels.append([p, y * a[1], Rad * a[2], A* a[3]])
                    with merged_count.get_lock():
                        merged_count += 1
                    del partials[f[1]]
                else:
                    partials[f[1]] = [f[0], y, Rad, A]
        if i % m == 0:
            lRels = len(Rels)
            lt = time.time()
            td = lt - ltd
            ltd = lt
            eta = td * (((ld / m) + (required_relations / lRels)) / 2)
            tds = humanfriendly.format_timespan(td)
            etas = humanfriendly.format_timespan(eta)
            msg = "relations_find: range(%d, %d), inverval: %d of %d, found: %d of %d, merged: %d, iter_elapsed: %s, eta: %s.\n" % (start,stop,i,(stop-start),lRels,required_relations, merged_count,tds,etas)
            sys.stderr.write(msg)
    
    td = time.time() - st
    tds = humanfriendly.format_timespan(td)
    msg = "relations_find: Ended range(%d, %d), inverval: %d of %d, found: %d of %d, merged: %d, time elapsed: %s\n" % (start,stop,i,(stop-start),len(Rels),required_relations, merged_count, tds)
    sys.stderr.write(msg)
    
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
    return marks, A

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
    primes, mod_root, log_p, num_prime, p = [], [], [], 0, 3
    while p < bound or num_prime < 3:
        leg = legendre(n % p, p)
        if leg == 1:
            primes += [p]
            #mod_root += mod_sqrt(n ,p)
            log_p += [log10(p)]
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
    sys.stderr.write("min_prime: %d, fudge: %f\n" % (min_prime,fudge))
    thresh -= fudge
    return min_prime, thresh, fudge


def generate_polys(N, Prime_base, x_max, needed):
    """ 
    It searchs for distinct needed polys congruent to n. 
    """
    n=1
    cpolys = 0
    polys = []
    polys_ABC = []
    early_factors = []
    while cpolys <= needed:
        pol = Poly(N, Prime_base, x_max, search = n, verbose=False)
        #print(pol.early_factors)
        if len(pol.early_factors) > 0:
            for early_factor in pol.early_factors:
                if early_factor not in early_factors:
                    early_factors.append(early_factor)
        else:
            pol_ABC = (pol.A,pol.B,pol.C)
            if pol_ABC not in polys_ABC:
                m = "Found poly: f(x) = %dx^2 + %dx + %d\n" % (pol_ABC)
                sys.stderr.write(m.replace("+ -","- "))
                polys_ABC.append(pol_ABC)
                polys.append(pol)
                cpolys += 1
        n += 1
    return polys, early_factors


def _MPQS(N, verbose=True, M = 1):
    """ 
    Main MPQS function. 
    """
    bN, lN = int(log2(N)), len(str(N))
    i2N = isqrt(N) 
    i2Np1 = i2N + 1 
    root_2n = isqrt(2*N)
    Rels = []

    T = cpu_count()

    B1 = pow(int(log10(pow(N, 6))),2) 
    Prime_base, log_p = find_primebase(N, B1)
    B1 //= M
    B2 = len(Prime_base)
  
    if B2 == 1:
        sys.stderr.write("Found small factor: %d\n" % Prime_base[0])
        return Prime_base + _MPQS(N // Prime_base[0])

    multiplier = choose_multiplier(N, Prime_base)
    Nm = multiplier * N
    
    x_max = B2 *60  # size of the sieve
    m_val = (x_max * root_2n) >> 1
    thresh = log10(m_val) * 0.735
    min_prime = int(thresh * 3)
    
    required_relations_ratio = 1.05
    required_relations = round(len(Prime_base) * required_relations_ratio) 

    sys.stderr.write("Factoring N: %d, bits: %d, digits: %d, B1: %d, B2: %d\n" % (N,bN,lN,B1,B2))
    sys.stderr.write("Multiplier is: %d\n" % multiplier)
    sys.stderr.write("Need %d relations\n" % (required_relations))

    start = 0
    stop = B1 # range to sieve

    polys, early_factors = generate_polys(Nm, Prime_base, x_max, T) # generate n distinct polys one for each cpu core.
    if len(early_factors) > 0:
        tmp = 1
        small = []
        for early_factor in early_factors:
            g = gcd(early_factor, N)    
            if N > g > 1: 
                tmp *= g
                small.append(g)
        if tmp > 1:
            return small + MPQS(N // tmp)
    
    manager = Manager()
    Rels = manager.list() # placeholder list for relations shareable between child processes.

    merged_count = manager.Value("i", 0)

    while True:
        # trim primes, recalc min
        Prime_base = [p for p in Prime_base if p > min_prime]
        smooth_base = prod(Prime_base)
        min_prime, thresh, fudge = recalculate_min_prime_thresh(thresh, Prime_base, log_p)

        t1 = time.time()
        sys.stderr.write("Data collection with %d threads...\n" % T)

        inputs = []
        
        # generate tasks parameters
        for poly in polys:
            s1 = min(poly.start_vals[0]) 
            s2 = max(poly.start_vals[0]) 
            inputs += [(Nm, start + s1, stop + s1 , Prime_base, smooth_base, Rels, merged_count, required_relations,  poly)]

        # deploy tasks to every cpu core.
        pols = []
        workpool = Pool(T)
        with workpool:
            R = workpool.starmap(relations_find, inputs)  
        workpool.close()
        workpool.join() 
        
        t2 = time.time()
        sys.stderr.write("Done in: %f secs.\n" % (t2-t1))
        sys.stderr.write("Found %d rels with %d base primes.\n" % (len(Rels),len(Prime_base)))

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
            
            if result != None:
                return result

        start = stop
        stop *= 2

        sys.stderr.write("Need to sieve %d differences more...\nNew sieving range %d:%d\n" % (required_relations - len(Rels),start,stop))
      

def MPQS(N):
    """ 
    Iterative version of MPQS. 
    """
    result = _MPQS(N)
    if result != None:
        R= []
        for r in result:
            if is_prime(r):
                sys.stderr.write("Found factor: %d\n" % r)
                R.append(r)
            else:
                R += _MPQS(r)
        return R


if __name__ == "__main__":
    N = Integer(sys.argv[1])
    ts = time.time()
    r = MPQS(N)
    print(r)
    td = time.time() - ts
    tds = humanfriendly.format_timespan(td) 
    sys.stderr.write("All done in: %s.\n" % (tds))
