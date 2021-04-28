#!/usr/bin/env python3
# Author Dario Clavijo 2021
# adaptation from https://en.wikipedia.org/wiki/Shor%27s_algorithm

import random
import sys
from gmpy2 import gcd, isqrt,log2

# This function is actualy not being called, its shown just for documenting purposes.
def multiplicativeOrder(a,n):
   """ Naive multiplicative_order implementation of orders of magnitude slower than sage's
   Algo complexity is O(N)"""
   
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
  """ Shor's algorithm recreated in a classical computer
  Premise: given a ^ P - 1 % N  == 0, where a is random, and p = multOrd(a,N) we get a factor o N.
  Sage's implementation of multiplicative order first factors the integer (N) in question with PARI.
  In other words it can be said that the complexity of multiplicativeOrder problem 
  is equal to the complexity of the factorization problem.
  This makes this algorithm only a showcase for studying purposes and not an improvement of performance in any way."""

  n = 1
  p = q = None
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
        if explain:
          print("pow(%d, %d ,%d) = %d" % (a,r//2,N,ar2)) 
        if N-1 > ar2 > 1 and (ar2 + 1) % N != 0:
          p,q = gcd(int(ar2+1),N),gcd(int(ar2-1),N)
          #print(x,y)
          if N > p > 1:
            if explain:
              print("gcd(%d ^ %d - 1, %d) = %d" % (a,r//2,N,p))
            q = N//p
          if N > q > 1:
            if explain:
              print("gcd(%d ^ %d - 1, %d) = %d" % (a,r//2,N,q))
            x = N//q
          if p != None and q != None:
            assert p*q == N
            return p,q
    n+=1

if __name__ == "__main__":
  """
  time sage shors_classic.sage 6611945602629820434117612619247836840607270940689358040985518679439408007369866558963727598868092310222457775451141333970448929023239
  log2(N) = 442, digits(N) = 133
  Iter: 1
  gcd(648364926028922147965451188485987542837109702292308105142056428695249584100146182958584631176408868457221639858637637175557795606755, 6611945602629820434117612619247836840607270940689358040985518679439408007369866558963727598868092310222457775451141333970448929023239) = 1
  multOrd(648364926028922147965451188485987542837109702292308105142056428695249584100146182958584631176408868457221639858637637175557795606755, 6611945602629820434117612619247836840607270940689358040985518679439408007369866558963727598868092310222457775451141333970448929023239) = 3379620273251187101376818688869414473378495161760299772869915717699712138568616296446109283202225832492363688900385902720
  pow(648364926028922147965451188485987542837109702292308105142056428695249584100146182958584631176408868457221639858637637175557795606755, 1689810136625593550688409344434707236689247580880149886434957858849856069284308148223054641601112916246181844450192951360 ,6611945602629820434117612619247836840607270940689358040985518679439408007369866558963727598868092310222457775451141333970448929023239) = 822680968913290200259975302793050771720099056310233973126016591324403552685984372332305899290644971490603783815121396422325799529066
  gcd(648364926028922147965451188485987542837109702292308105142056428695249584100146182958584631176408868457221639858637637175557795606755 ^ 1689810136625593550688409344434707236689247580880149886434957858849856069284308148223054641601112916246181844450192951360 - 1, 6611945602629820434117612619247836840607270940689358040985518679439408007369866558963727598868092310222457775451141333970448929023239) = 93461639715357977769163558199606896584051237541638188580280321
  N = 70745020339540650970024908966159700863734087283476256855941113111428359 * 93461639715357977769163558199606896584051237541638188580280321

  real	0m3.648s
  user	0m3.571s
  sys	0m1.054s
  """
  N=Integer(sys.argv[1])
  if len(sys.argv) > 2:
    scaler = 10**int(sys.argv[2])
  print("log2(N) = %d, digits(N) = %d" % (N.nbits(),N.ndigits()))
  p,q=shor_factor(N,explain=True) 
  print("N = %d * %d" % (p,q))     
