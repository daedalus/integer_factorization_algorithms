#!/usr/bin/env python3
# Multiple polynomial Quadratic sieve
# Author Dario Clavijo 2021
# License: GPLv3

import time
import sys
from gmpy2 import gcd, gcdext, isqrt, is_prime, next_prime, log2, log10, legendre, powmod, invert
from sage.parallel.multiprocessing_sage import parallel_iter
from multiprocessing import cpu_count, Pool, Manager
from itertools import repeat
import humanfriendly

def mod_sqrt(n, p):
    """ tonelli shanks algorithm """
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


def filter_p(ppws):
    """ filter out even powers """
    d = {}
    for p,pw in ppws:
        if p not in d:
            d[p] = pw
        else:
            d[p] += pw
    #print(d)
    tmp2 = []
    for p in d:
        if d[p] & 1 != 0: # only keep odd powers
            tmp2.append(p)
    return tmp2


def trial_division_minus_even_powers(n, P):
    a, r, n = trial_division(n, P)
    a = filter_p(a)
    return [a,r,n]


def is_smooth(x, P):
    """ check if n is B-smooth """
    y = x
    for p in P:
        while y % p ==0:
            y //= p
    return abs(y) == 1 


def minifactor2(x, P, D = 1):
    """ tries to factor out common primes with gcd and filter out even power primes"""
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
        tmp = filter_p(a)  
        return (tmp, R, x), D                 


def minifactor3(x, P, smooth_base):
    """ minifactor """
    x = abs(x)
    smooth = gcd(x, smooth_base)
    if smooth > 1:
        not_smooth = x // smooth
        a1, r1, n1 = trial_division(smooth, P)
        a2, r2, n2 = trial_division(not_smooth, P)
        if r2 == 1:
            a = a1 + a2
            a3 = filter_p(a)
            return a3, r2, x
        
                         
def minifactor(x, P):
    """ minifactor algo, finds odd-power primes in composites """
    p = trial_division_minus_even_powers(x,P)
    if p[1] == 1: # if 1 x is B-smooth
        return p


class Poly:
    """ Quadratic polynomial helper class: 
    type Ax + 2Bx + C """
    def __init__(self, n, P, x_max, search = 0, verbose = None):
        self.n = int(n)
        self.P = P
        self.x_max = x_max
        self.search = search
        self.verbose = verbose
        self.create()
        self.solve_for_x()
 
    def create(self):
        n = self.n
        x_max = self.x_max
        root_2n = isqrt(2*n)
        root_A = next_prime(isqrt(root_2n // x_max))
        self.root_A = root_A
        s=0
        while True:
            root_A = next_prime(root_A)
            leg = legendre(n, root_A)
            if leg == 1:
                if s > self.search:
                    break
            elif leg == 0:
                self.early_factor = root_A
                self.A = None
                self.B = None
                self.C = None
                return root_A, None, None
            s+=1

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
        A = self.A
        B = self.B
        C = self.C
        Ax = A * x
        Rad = Ax + B 
        return (Rad + B) * x + C, Rad # same as (Ax + 2B)x + C
        


def relations_find(N, start, stop, P, smooth_base, Rels, required_relations, pol = None):
    """ relations search funcion """
    #print(N,start,stop)
    sys.stderr.write("relations_find: range(%d, %d), interval: %d sieving start\n" % (start, stop, (stop-start)))
    if (stop-start) < 0:
        return [] 

    I = [x for x in range(start, stop) if not is_prime(x) and not is_square(x)]
   
    #D = reduce(lambda x, y: x * y, I)
    Diffs = [(pol.eval(x),x) for x in I]
    Found_Rels = []

    A = pol.A
    #B = pol.B
    #C = pol.C

    ltd = time.time()
    ld = len(Diffs)
    m = 1000
    msg = ""
    for i in range(ld):
        if len(Rels) > required_relations:
            break
        yRad, x = Diffs[i]
        y, Rad = yRad
        #r = minifactor2(y, P, D//y)
        #f = minifactor(y, P)
        f = minifactor3(y, P, smooth_base)

        #print(r,trial_division(y,P),y)
        if 1:
        #if r != None:
            #f, D = r
            if f != None and f[1] == 1:
                Rels.append((f,(y,Rad,A,x)))
        if i % m == 0:
            lt = time.time()
            td = lt - ltd
            ltd = lt
            eta = td * (ld/m)
            tds = humanfriendly.format_timespan(td)
            etas = humanfriendly.format_timespan(eta)
            msg = "relations_find: range(%d, %d), inverval: %d of %d, found: %d of %d, iter_elapsed: %s, eta: %s.\n" % (start,stop,i,(stop-start),len(Rels),required_relations,tds,etas)
            sys.stderr.write(msg)
    sys.stderr.write(msg)
     
    #Rels += Found_Rels
    
    return Found_Rels


def linear_algebra(Rels, P):
    """ linear algebra, it generates a matrix in GF(2) then computes it's null-space"""
    M = matrix(GF(2), len(Rels), len(P), lambda i, j:P[j] in Rels[i][0][0])
    return M.left_kernel().basis()


def process_basis_vectors(N, basis, Rels):
    """ process each basis vector """
    for K in basis:
        lhs = rhs = Ahs = 1
        I = [f for f, k in zip(Rels, K) if k == 1]
        if len(I) > 0:
            for i in I:
                lhs *= i[1][0] # left-hand side
                rhs *= i[1][1] # rigth-hand side
                Ahs *= i[1][2] # A-term in poly
            LHS = Ahs * lhs
            g = gcd(isqrt(LHS)-rhs,N)
            if N > g > 1:
                return [int(g), int(N//g)]


def find_primebase(n, bound):
    """ finds the base prime for given n and bound """
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
    """ recalculates the min efective prime base bound """
    min_prime = int(thresh * 3)
    fudge = sum(log_p[i] for i, p in enumerate(Prime_base) if p < min_prime)
    fudge = fudge // 4
    sys.stderr.write("min_prime: %d, fudge: %f\n" % (min_prime,fudge))
    thresh -= fudge
    return min_prime, thresh, fudge


def generate_polys(N, Prime_base, x_max, needed):
    """ It searchs for distinct needed polys congruent to n."""
    n=1
    cpolys = 0
    polys = []
    polys_ABC = []
    early_factor = None
    while cpolys <= needed:
        pol = Poly(N, Prime_base, x_max, search=n, verbose=False)
        if pol.A == None and pol.B == None and pol.C == None:
            early_factor = pol.early_factor
            polys = None
            break
        n += 1
        pol_ABC = (pol.A,pol.B,pol.C)
        if pol_ABC not in polys_ABC:
            m = "Found poly: f(x) = %dx^2 + %dx + %d\n" % (pol_ABC)
            sys.stderr.write(m.replace("+ -","- "))
            polys_ABC.append(pol_ABC)
            polys.append(pol)
            cpolys += 1
    return polys, early_factor


def _MPQS(N, verbose=True, M = 1):
    """ main MPQS function """
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

    x_max = B2 *60  # size of the sieve

    m_val = (x_max * root_2n) >> 1
    thresh = log10(m_val) * 0.735
    min_prime = int(thresh * 3)
    
    required_relations_ratio = 1.05
    required_relations = round(len(Prime_base) * required_relations_ratio) 

    sys.stderr.write("Factoring N: %d, bits: %d, digits: %d, B1: %d, B2: %d\n" % (N,bN,lN,B1,B2))
    sys.stderr.write("Need %d relations\n" % (required_relations))

    start = 0
    stop = B1

    polys, early_factor = generate_polys(N, Prime_base, x_max, T) # generate n distinct polys one for each cpu core.
    if polys == None and early_factor != None:
        sys.stderr.write("Found small factor: %d\n" % early_factor)
        return [early_factor] + _MPQS(N // early_factor)

    manager = Manager()
    Rels = manager.list()

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
            inputs += [(N, start + s1, stop + s1 , Prime_base, smooth_base, Rels, required_relations,  poly)]

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
            
            basis = linear_algebra(Rels, Prime_base)

            t3 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t3-t2))
            sys.stderr.write("Matrix reduction...\n")
            result = process_basis_vectors(N,basis, Rels)
            t4 = time.time()
            sys.stderr.write("Done in: %f secs.\n" % (t4-t3))
            
            if result != None:
                return result

        start = stop
        stop *= 2

        sys.stderr.write("Need to sieve %d differences more...\nNew sieving range %d:%d\n" % (required_relations - len(Rels),start,stop))
      

def MPQS(N):
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
    td = time.time()
    sys.stderr.write("All done in: %f secs.\n" % (td-ts))
