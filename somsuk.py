#!/usr/bin/env python3
# Author Dario Clavijo 2022

from gmpy2 import isqrt, gcd, sqrt, powmod, isqrt_rem, is_square
import random
import time

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


def somsuk(n):
  """
  Implementation based on Kritsanapong Somsuk paper: 
  The new integer factorization algorithm based on fermat’s
  factorization algorithm and euler’s theorem (2020).
  ISSN: 2088-8708, DOI: 10.11591/ijece.v10i2.pp1469-1476
  """
  u = isqrt(n) << 1
  c = n
  while gcd(c, n) > 1:
    c = random.randrange(2, n - 1)
  a = powmod(c, -1, n)
  s = powmod(c, 2, n)
  phi_aprox = n - u  + 1
  t = powmod(a, phi_aprox, n)
  while True:
    if t == 1:
      x = u >> 1
      y = isqrt(x * x - n)
      return x - y, x + y
    t = (t*s) % n
    u += 2

def timeit(f,n):
  t0 = time.time()
  r = f(n)
  t1 = time.time()
  print(str(f).split(" ")[1],"(", n , ") =", r,"time:",round(t1-t0,5))

def test():
  timeit(somsuk, 340213)
  timeit(fermat, 340213)
  timeit(somsuk, 788582867650121563)
  timeit(fermat, 788582867650121563)
  timeit(somsuk, 1047329636821139813)
  timeit(fermat, 1047329636821139813)

if __name__ == "__main__":
  test()


