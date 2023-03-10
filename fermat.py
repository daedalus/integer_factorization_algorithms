import sys
from gmpy2 import isqrt, isqrt_rem, is_square

def fermat(n):
    a, rem = isqrt_rem(n)
    b2 = -rem
    c0 = (a << 1) + 1
    c = c0
    while not is_square(b2):
        b2 += c
        c += 2
    a = (c - 1) >> 1
    b = isqrt(b2)
    return a - b, a + b

print(fermat(int(sys.argv[1])))
