from gmpy2 import *

def factor(N):
  """
  Algorithm described in the wagstaf-joy of factoring book.
  """
  f = N1 = N - 1
  while f & 1 == 0: f >>= 1
  a = 2
  while a <= N1:
    r1 = powmod(a, f << 1, N) 
    if r1 == 1:
      r = powmod(a, f, N) 
      p = gcd(r - 1, N)
      q = gcd(r + 1, N)
      if (q > p > 1): #and (p * q == N):
        return factor(p) + factor(q)
    a = next_prime(a)
  return [N]

def tests():
  for n in [561,6601, 8911,23224518901,1901 * 3701,6601 * 8911,52193 * 3701]:
    print(n, factor(n))

if __name__ == "__main__":
  tests()
