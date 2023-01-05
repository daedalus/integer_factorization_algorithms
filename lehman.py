from gmpy2 import *


def cuberoot(n):
  a,b = iroot(n,3)
  return a


def lehman(n):
    """
    based on: https://programmingpraxis.com/2017/08/22/lehmans-factoring-algorithm/
    """
    if is_congruent(n, 2, 4):
        return []

    for k in range(1, cuberoot(n), 1):
        nk4 = n*k << 2
        ki4 = isqrt(k) << 2
        ink4 = isqrt(nk4) + 1
        i6,_ = iroot(n, 6)
        ink4i6ki4 = ink4 + (i6 // (ki4))
        for a in range(ink4, ink4i6ki4 + 1, 2):
            b2 = (a * a) - nk4
            if is_square(b2):
                b = isqrt(b2)
                p = gcd(a + b, n)
                q = gcd(a - b, n)
                return p, q
    return []

def tests():
  N = [29 * 89, 3*11, 3141592651*3141592661, 7919*10861]
  N += [94738740796943840961823530695778701408987757287583492665919730017973847138345511139064596113422435977583856843887008168202003855906708039013487349390571801141407245039011598810542232029634564848797998534872251549660454277336502838185642937637576121533945369150901808833844341421315429263207612372324026271327]
  N += [9733382803370256893136109840971590971460094779242334919432347801491641617443615856221168611138933576118196795282443503609663168324106758595642231987246769*9733382803370256893136109840971590971460094779242334919432347801491641617443615856221168611138933576118196795282443503609663168324106758595642231987248907]
  N=sorted(N)
  for n in N:
    p,q = lehman(n)
    print(n,p,q)
    assert p*q == n

if __name__ == "__main__":
  tests()

