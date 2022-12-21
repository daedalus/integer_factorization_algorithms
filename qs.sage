#!/usr/bin/env python3
# Code taken from https://facthacks.cr.yp.to/facthacks.sage
# FactHacks is joint work by Daniel J. Bernstein, Nadia Heninger, and Tanja Lange.

from gmpy2 import gcd, isqrt, is_prime, next_prime
import itertools


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


def minifactor(x, P):
    if not is_prime(x) and not is_square(x):
        p = trial_division_minus_even_powers(x, P)
        if p[1] == 1: # if 1 x is B-smooth
            return p


def QS(N,differences,y):
    i2N = isqrt(N)
    X = [int(a * a - N) for a in range(i2N + 1, i2N + differences)]
    P = list(primes(2, y))
    F = list(filter(None, map(minifactor, X, itertools.repeat(P, len(X)))))

    M = matrix(GF(2), len(F), len(P), lambda i, j:P[j] in F[i][0])
    for K in M.left_kernel().basis():
        I = [f for f, k in zip(F, K) if k == 1]
        x = product([int(isqrt(f[2] + N)) for f in I])
        y = isqrt(product([f[2] for f in I]))
        g0, g1 = gcd(N, x - y), gcd(N, x + y)
        if N > g0 > 1 or N > g1 > 1:
            return g0, g1


if __name__ == "__main__":
    N = int(sys.argv[1])
    y ,_ ,diff = prebuilt_params(N.bit_length())
    r = QS(N, diff * 10, y * 10)
    if r != None:
        print("%d, %d" % r)

