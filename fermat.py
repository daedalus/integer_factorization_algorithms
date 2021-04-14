import sys
from gmpy2 import isqrt

def fermat(n):
        """Fermat attack"""
        a = isqrt(n)
        a = b
        b2 = (a ** 2) - n
        count = 0
        while (b ** 2) != b2:
            a = a + 1
            b2 = (a ** 2) - n
            b = isqrt(b2)
            count += 1
        p = a + b
        q = a - b
        assert n == p * q
        return p, q

print(fermat(int(sys.argv[1])))
