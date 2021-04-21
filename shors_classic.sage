#!/usr/bin/env python3
# Author Dario Clavijo 2021
# adaptation from https://en.wikipedia.org/wiki/Shor%27s_algorithm

import random
import sys
from gmpy2 import gcd, isqrt

# This function is actualy not being called, its shown just for documenting purposes.
def multiplicativeOrder(a,n):
   """ Naive multiplicative_order implementation of orders of magnitude slower than sage's """
   """ Algo complexity is O(N) """
   if gcd(a,n) > 1:
     return -1
   k = 2
   while k < n:
     p = pow(a,k,n)
     if p == 1:
       return k
     k += 1
   return -1

def shor_factor(N,explain=False):
  """ Shor's algorithm recreated in a classical computer """
  """ a ^ P - 1 % N  == 0, where a is random, and p = multOrd(a,N) """
  n=1
  while True:
    a = random.randint(2,N-1)
    b = gcd(a,N)
    if explain:
      print("Iter: %d" % n)
      print("gcd(%d, %d) = %d" % (a,N,b))
    if 1 < b < N:
      return N//b,b
    else:
      r = Mod(a,N).multiplicative_order() # Equivalent to quantum part of the algorithm, complexity: O(N).
      if explain:
        print("multOrd(%d, %d) = %d" % (a,N,r))
      if r % 2 == 0:
        ar2 = pow(a,(r//2),N)  
        if (ar2 + 1) % N != 0:
          x,y = gcd(int(ar2+1),N),gcd(int(ar2-1),N)
          #print(x,y)
          if N > x > 1:
            if explain:
              print("gcd(%d ^ %d - 1, %d) = %d" % (a,r//2,N,x))
            return int(N//x),int(x)
          if N > y > 1:
            if explain:
              print("gcd(%d ^ %d - 1, %d) = %d" % (a,r//2,N,y))
            return int(N//y),int(y)
        else:
            print("nok")
    n+=1

N=int(sys.argv[1])
p,q=shor_factor(N,explain=True) 
print("N = %d * %d" % (p,q))     
