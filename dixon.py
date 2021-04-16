from gmpy2 import isqrt,gcd,next_prime
import sys

def dixon(N,B=7):

  def primes(B):
    p = 2
    tmp = [p]
    while p < B-1:
      p = next_prime(p)
      tmp.append(p)
    return tmp

  base = primes(B)
  start = isqrt(N)
  i = start

  basej2N = []
  for j in range(0,len(base)):
    basej2N.append(pow(base[j],2,N))

  while i < N:
    i2N = pow(i,2,N)
    for k in range(0,len(base)):
      if i2N == basej2N[k]:
        f=gcd(i - base[k],N)
        if 1 < f < N:
          return f,N//f
    i+=1
  return -1

print(dixon(int(sys.argv[1])))
