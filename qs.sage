#!/usr/bin/env python3
# Code taken from https://facthacks.cr.yp.to/facthacks.sage
# FactHacks is joint work by Daniel J. Bernstein, Nadia Heninger, and Tanja Lange.

from gmpy2 import gcd, isqrt, is_prime, next_prime
import itertools


def trial_division_minus_even_powers(n,P):
    a = []
    r = n
    pw = 0
    if r & 1 == 0:
        while r & 1 == 0:
            pw += 1
            r >>= 1
        if pw & 1 != 0:
            if 2 not in a:
                a.append(2)
    l=len(P)
    i=0
    while (i < l) and (P[i] <= n):
        if r % P[i] == 0:
            pw = 0
            while r % P[i] == 0:
                pw += 1
                r //= P[i]
            if pw & 1 != 0: 
                if P[i] not in a:
                    a.append(P[i])
        else:
            i += 1
    return [a,r,n]


def minifactor(x,P):
    if not is_prime(x):
        i = isqrt(x)
        if pow(i,2) != x: # perfect squares are even powers
            p = trial_division_minus_even_powers(x,P)
            if p[1] == 1: # if 1 x is B-smooth
                return p


def QS(N,differences,y):
    i2N = isqrt(N)
    X = [Integer(a^2 - N) for a in range(i2N + 1,i2N + differences)]
    P = list(primes(2, y))
    F = list(filter(None,map(minifactor, X, itertools.repeat(P, len(X)))))

    M = matrix(GF(2), len(F), len(P), lambda i, j:P[j] in F[i][0])
    for K in M.left_kernel().basis():
        I = [f for f, k in zip(F, K) if k == 1]
        x = product([int(isqrt(f[2] + N)) for f in I])
        y = isqrt(product([f[2] for f in I]))
        g0, g1 = gcd(N, x - y), gcd(N, x + y)
        if N > g0 > 1:
            return g0, N // g0
        if N > g1 > 1:
            return g1, N // g1


if __name__ == "__main__":
    r = QS(Integer(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
    if r != None:
        print("%d, %d" % r)

