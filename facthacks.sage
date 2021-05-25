# code taken from: https://facthacks.cr.yp.to/facthacks.sage

def greatestcommondivisor(x,y):
  while x != 0: x,y = y%x,x
  return abs(y)

# example:
print greatestcommondivisor(4187,5989)
# output: 53
# double check:
print gcd(4187,5989)
# output: 53

def product(X):
  if len(X) == 0: return 1
  while len(X) > 1:
    X = [prod(X[i*2:(i+1)*2]) for i in range((len(X)+1)/2)]
  return X[0]

# example:
print product([314,159,265,359,897])
# output: 4260489878970
# double check:
print 314*159*265*359*897
# output: 4260489878970

def producttree(X):
  result = [X]
  while len(X) > 1:
    X = [prod(X[i*2:(i+1)*2]) for i in range((len(X)+1)/2)]
    result.append(X)
  return result

# example:
print producttree([10,20,30,40,50,60])
# output: [[10, 20, 30, 40, 50, 60], [200, 1200, 3000], [240000, 3000], [720000000]]

def remaindersusingproducttree(n,T):
  result = [n]
  for t in reversed(T):
    result = [result[floor(i/2)] % t[i] for i in range(len(t))]
  return result

def remainders(n,X):
  return remaindersusingproducttree(n,producttree(X))

# example:
print remainders(8675309,[11,13,17,19,23])
# output: [5, 6, 5, 4, 8]
# double check:
print [8675309 % p for p in [11,13,17,19,23]]
# output: [5, 6, 5, 4, 8]

def batchgcd_simple(X):
  R = remainders(product(X),[n^2 for n in X])
  return [gcd(r/n,n) for r,n in zip(R,X)]

def batchgcd_faster(X):
  prods = producttree(X)
  R = prods.pop()
  while prods:
    X = prods.pop()
    R = [R[floor(i/2)] % X[i]**2 for i in range(len(X))]
  return [gcd(r/n,n) for r,n in zip(R,X)]

# example:
print batchgcd_simple([1909,2923,291,205,989,62,451,1943,1079,2419])
print batchgcd_faster([1909,2923,291,205,989,62,451,1943,1079,2419])

def primesin(P,x):
  result = remainders(x,P)
  return [p for p,r in zip(P,result) if r == 0]

def primesinproduct(P,X):
  return primesin(P,product(X))

def primesineach(P,X):
  n = len(X)
  if n == 0: return []
  P = primesinproduct(P,X)
  if n == 1: return [P]
  return primesineach(P,X[:n//2]) + primesineach(P,X[n//2:])

# example:
P = list(primes(2,10))
X = [50,157,266,377,490,605]
print primesineach(P,X)
# output: [[2, 5], [], [2, 7], [], [2, 5, 7], [5]]
# double check:
print [[p for p in P if x%p==0] for x in X]
# output: [[2, 5], [], [2, 7], [], [2, 5, 7], [5]]

def fermatfactor(N):
  if N <= 0: return [N]
  if is_even(N): return [2,N/2]
  a = ceil(sqrt(N))
  while not is_square(a^2-N):
    a = a + 1
  b = sqrt(a^2-N)
  return [a - b,a + b]

print fermatfactor(91)
print fermatfactor(1009)
print fermatfactor(6557)
N = 115792089237316195423570985008721211221144628262713908746538761285902758367353
print fermatfactor(N)
N = 115792089237316195448679392282006640413199890130332179010243714077028592474181
print fermatfactor(N)

def lattice(n,nearp,howclose,t,k):
  R.<x> = PolynomialRing(ZZ)
  f = howclose*x+nearp
  M = matrix(t)
  for i in range(t):
    M[i] = (f^i*n^max(k-i,0)).coeffs()+[0]*(t-1-i)
  M = M.LLL()
  Q = sum(z*(x/howclose)^i for i,z in enumerate(M[0]))
  for r,multiplicty in Q.roots():
    if nearp+r > 0:
      g = gcd(n,nearp+r)
      if g > 1: return [g,n/g]
  return [1,n]

# examples:
print lattice(314159265358979323,317213000,1000,5,2)
# output: [317213509, 990371647]
print lattice(314159265358979323,317210000,10000,40,20)
# output: [317213509, 990371647]

def oddpowers(P,x):
  if len(P) == 0: return [[],x,x]
  P = primesin(P,x)
  e,r,x = oddpowers([p*p for p in P],x)
  Q = primesin(P,r)
  return [Q,r / product(Q),x]

def oddpowersineach(P,X):
  return [oddpowers(Q,x) for Q,x in zip(primesineach(P,X),X)]

def easilyfactorable(x,pmodx):
  smoothpart = Mod(pmodx,x)^x.nbits()
  return (x / gcd(x,Integer(smoothpart))).is_square()

def easyfactorizations(P,X):
  result = zip(X,remainders(product(P),X))
  result = [x for x,pmodx in result if easilyfactorable(x,pmodx)]
  return oddpowersineach(P,result)

# example:
print easyfactorizations([2,3,5,7],[50,157,266,377,490,605])
# output: [[[2], 1, 50], [[2, 5], 1, 490], [[5], 121, 605]]

def qs_basic(N,differences,y):
  X = [Integer(a^2-N) for a in range(sqrt(N)+1,sqrt(N)+differences)]
  P = list(primes(2,y))
  F = easyfactorizations(P,X)
  M = matrix(GF(2),len(F),len(P),lambda i,j:P[j] in F[i][0])
  for K in M.left_kernel().basis():
    x = product([sqrt(f[2]+N) for f,k in zip(F,K) if k==1])
    y = sqrt(product([f[2] for f,k in zip(F,K) if k==1]))
    return [gcd(N,x - y),gcd(N,x + y)]
  return [1,N]

# example:
print qs_basic(275801,1000,20)
# output: [389, 709]

def rho(n):
  # Pollard's rho method
  c = 10
  a0 = 1
  a1 = a0^2+c
  a2 = a1^2+c
  while gcd(n,a2-a1) == 1:
    a1 = (a1^2+c) % n
    a2 = (a2^2+c) % n
    a2 = (a2^2+c) % n
  g = gcd(n,a2-a1)
  return [g,n / g]

# examples (second example is faster because prime is smaller):
print rho(314159265358979323)
# output: [990371647, 317213509]
print rho(698599699288686665490308069057420138223871)
# output: [2053, 340282366920938463463374607431768211507]

def p1exponent(cutoff):
  return lcm(range(1,cutoff))

def p1(n,cutoff):
  # Pollard's p-1 method
  g = gcd(n,Integer(pow(2,p1exponent(cutoff),n))-1)
  return [g,n/g]

# example:
print p1(38568900844635025971879799293495379321,2^14)
