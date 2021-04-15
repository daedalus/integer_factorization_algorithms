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

  i2N = []
  while i<=N:
    i2N.append(pow(i,2,N))
    i+=1

  basej2N = []
  for j in range(0,len(base)):
    basej2N.append(pow(base[j],2,N))

  for i in range(0,len(i2N)):
    for k in range(0,len(base)):
      if i2N[i] == basej2N[k]:
        f=gcd(start + i - base[k],N)
        if 1 < f < N:
          return f,N//f
  return -1

print(dixon(int(sys.argv[1])))
