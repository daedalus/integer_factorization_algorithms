from gmpy2 import isqrt, gcd, is_prime
import sys

def pollard_rho(self, n, seed=2, p=2, c=1):
        if is_prime(n):
            return n
        if n % 2 == 0:
            return 2
        if n % 3 == 0:
            return 3
        if n % 5 == 0:
            return 5
        f = lambda x: x ** p + c
        x, y, d = seed, seed, 1
        while d == 1:
            x = f(x) % n
            y = f(f(y)) % n
            d = gcd((x - y) % n, n)
            if N > d > 1:
                return d
        
print(int(sys.argv[1]))
