from sympy.core.power import isqrt
from sympy import isprime
from functools import cache

@cache
def factorMersenne(n):
  " Lazzy programming implementation "
  if isprime(n):  yield n
  p = n.bit_length()
  for q in range(2*p + 1, isqrt(n) + 1, 2*p):
    if [0,1,0,0,0,0,0,1][q % 8] and pow(2, p, q) == 1:
      if n % q == 0:
        yield from factorMersenne(q)
        yield from factorMersenne(n//q)

def g():
  for x in range(1, 939):
    n = (1 << x)-1
    print(x, n,[p for p in factorMersenne(n)])
    print(factorMersenne.cache_info())

g()
