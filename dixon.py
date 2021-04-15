from gmpy2 import isqrt,gcd,next_prime
import sys

def dixon(N,B=100):

  tmp = []
  pairs = []

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
  while i<=N:
    for j in range(len(base)):
      l = pow(i,2, N)
      r = pow(base[j],2,N)
      if l == r:
        pairs.append([i,base[j]])
        print(pairs)
    i+=1

  for i in range(len(pairs)):
    x = pairs[i][0]
    y = pairs[i][1]
    tmp.append(gcd(x-y,N))

  return tmp

print(dixon(int(sys.argv[1])))
