#!/usr/bin/env python3
# Code taken from https://facthacks.cr.yp.to/facthacks.sage
# FactHacks is joint work by Daniel J. Bernstein, Nadia Heninger, and Tanja Lange.

from gmpy2 import gcd, isqrt, is_prime, next_prime, is_power
import itertools
import sys


def trial_division_odd_pows(n, P):
    a, r, pw = [], n, 0
    while (r & 1) == 0:
        pw += 1
        r >>= 1
    if (pw & 1 == 1): a.append(2)
    for i in range(1, len(P)):
        pw, p = 0, P[i]
        while (r % p) == 0:
            pw += 1
            r //= p
        if (pw & 1 == 1): a.append(p)
    return a, r, n


def minifactor(x, P):
    p = trial_division_odd_pows(x, P)
    if p[1] == 1: return p # if 1 then x is b-smooth


def eulercriterion(N, p):
  return pow(N, (p - 1) >> 1, p)


def QS(N, B1):
    sys.stderr.write("max prime: %d\n" % B1)
    P = [p for p in [*(primes(2, B1))] if eulercriterion(N, p) == 1]
    i2n, B2, X, offset, lp = isqrt(N), 65536, [], 0, len(P)
    while True:
        S = i2n + offset
        sys.stderr.write("Relations matix: [cols: %d x rows: %d]\n" % (lp, offset + B2))
        X += [int(a * a - N) for a in range(S + 1, S + B2)]
        X = [*(filter(lambda x:not is_power(x) and not is_prime(x), X))]
        F = [*(filter(None, map(minifactor, X, itertools.repeat(P, len(X)))))]
        M = matrix(GF(2), len(F), lp, lambda i, j:P[j] in F[i][0])
        sys.stderr.write("performing linear algebra:\n")
        for K in M.left_kernel().basis():
            I = [f for f, k in zip(F, K) if k == 1]
            if len(I) > 0:
                sys.stderr.write("On Vector basis, taking square roots x and y.\n")
                x = product([isqrt(f[2] + N) for f in I])
                y = isqrt(product([f[2] for f in I]))
                g0, g1 = int(gcd(N, x - y)), int(gcd(N, x + y))
                if N > g0 > 1: return g0, N // g0
                if N > g1 > 1: return g1, N // g1
        offset += B2

if __name__ == "__main__":
    N = int(sys.argv[1])
    B1 = int(exp(sqrt(log(N) * log(log(N)))) ** (1 / sqrt(2)))  
    print(QS(N, B1))
