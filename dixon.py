from gmpy2 import isqrt,gcd,next_prime
import sys

def dixon(N,B=7):

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

  #print(base)
  
  start = isqrt(N)
  i = start
  while i<=N:
    for j in range(len(base)):
      l = pow(i,2, N)
      r = pow(base[j],2,N)
      #print(i,j,l,r)
      if l == r:
        pairs.append([i,base[j]])
        #print(pairs)
    i+=1

  #print(pairs)

  for i in range(len(pairs)):
    x = pairs[i][0]
    y = pairs[i][1]
    f = (gcd(x-y,N))
    if 1 < f < N:
      return f,N//f

print(dixon(int(sys.argv[1])))
