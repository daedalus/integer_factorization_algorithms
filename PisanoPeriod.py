#!/usr/bin/env python3
"""
Integer factorization with pisano period
Original repo https://github.com/wuliangshun/IntegerFactorizationWithPisanoPeriod/
White paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8901977
"""

import random
import time
from gmpy2 import *
import sys
sys.setrecursionlimit(5000)
from multiprocessing import Pool, cpu_count, Manager


class Fibonacci:
    def __init__(self, progress = False, verbose = True, multitasked = True):
        self.progress = progress
        self.verbose = verbose
        self.multitasked = multitasked

    def _fib_res(self,n,p):
        """ fibonacci sequence nth item modulo p """
        if n == 0:
            return (0, 1)
        else:
            a, b = self._fib_res(n >> 1,p)
            c = mod((mod(a, p) * mod(((b << 1) - a), p)), p)
            d = mod((powmod(a, 2, p) + powmod(b, 2, p)), p)
            if n & 1 == 0: 
                return (c, d)
            else:
                return (d, mod((c + d), p))


    def get_n_mod_d(self,n,d, use = 'mersenne'):
        if n < 0:
            ValueError("Negative arguments not implemented")
        if use == 'gmpy':
            return mod(fib(n), d)
        elif use == 'mersenne':
            return powmod(2, n, d) - 1
        else:
            return self._fib_res(n,d)[0]

    def ranged_period(self, N, start, stop, look_up_dest):
        print("ranged_period (%d,%d) start" % (start,stop))
        tmp_look_up = {}
        for x in range(start, stop):
            tmp_look_up[self.get_n_mod_d(x, N)] = x
        look_up_dest.update(tmp_look_up)
        #look_up_dest |= tmp_look_up
        print("ranged_period (%d,%d) end" % (start,stop))
        return 1
    

    def get_period_bigint(self, N, min_accept, xdiff, verbose = False):            
        search_len = int(pow(N, (1.0 / 6) / 100))
        
        if search_len < min_accept:
            search_len = min_accept
  
        if self.verbose:
            print('Search_len: %d, log2(N): %d' % (search_len,int(log2(N))))
        
        starttime = time.time()
        diff = xdiff 
        p_len = int((len(str(N)) + diff) >> 1) + 1
        begin = N - int('9'*p_len) 
        if begin <= 0:
            begin = 1
        end = N + int('9' * p_len)
    
        if self.verbose:    
            print('Search begin: %d, end: %d'%(begin, end))
        
        if self.multitasked and search_len > 1000:
            C = cpu_count() * 2
            search_work = search_len // C

            manager =  Manager()
            look_up = manager.dict()

            if self.verbose:
                print("Precompute LUT with %d tasks..." % C)

            inputs = []
            for x in range(0, search_len, search_work):
                inputs += [(N, x, x + search_work, look_up)]

            workpool = Pool(C)

            with workpool:
                results = workpool.starmap(self.ranged_period,inputs)
                print(results)
                workpool.close()
                workpool.join()

        else:
            look_up = {}
            self.ranged_period(N, 0, search_len, look_up)

        if self.verbose:
            print("LUT creation ended size: %d..." % len(look_up))
            print("Searching...")
        

        while True:       
            randi = random.randint(begin,end)            
            res = self.get_n_mod_d(randi, N)
            if res > 0:
                #print(res, res in look_up)
                if res in look_up:
                    res_n = look_up[res]
                    T = randi - res_n   
                    if T & 1 == 0: 
                        if self.get_n_mod_d(T, N) == 0:
                            td = int(time.time() - starttime)
                            if self.verbose:
                                print('For N = %d Found T:%d, randi: %d, time used %f secs.' % (N , T, randi, td))
                            return td, T, randi
                        else:
                            if self.verbose:
                                print('For N = %d\n Found res: %d, res_n: %d , T: %d\n but failed!' % (N, res, res_n, T))
            else:
                if randi & 1 == 0:
                    T = randi
                    td = int(time.time() - starttime)
                    if self.verbose:
                        print('First shot, For N = %d Found T:%d, randi: %d, time used %f secs.' % (N , T, randi, td))
                    return td, T, randi

  
    def _trivial_factorization_with_n_phi(self, N, T):
        phi = abs(N - T) + 1
        p1 = []
        d2 = N << 2

        phi2 = pow(phi, 2)
        phi2p4d = phi2 + d2
        phi2m4d = phi2 - d2

        if phi2m4d > 0:
            iphi2m4d = isqrt(phi2m4d)
            p1.append(phi + iphi2m4d)
            p1.append(phi - iphi2m4d)

        if phi2p4d > 0:
            iphi2p4d = isqrt(phi2p4d)
            p1.append(phi + iphi2p4d)
            p1.append(phi - iphi2p4d)

        if phi2m4d > 0:
            iphi2m4d = isqrt(phi2m4d)
            p1.append(-phi + iphi2m4d)
            p1.append(-phi - iphi2m4d)

        if phi2p4d > 0:
            iphi2p4d = isqrt(phi2p4d)
            p1.append(-phi + iphi2p4d)
            p1.append(-phi - iphi2p4d)

        for p in p1:
            g = gcd((p >> 1),N)
            if N > g > 1:
                return int(g),int(N//g)
   

    def factorization(self, N, min_accept, xdiff):
        res = self.get_period_bigint(N, min_accept, xdiff, verbose=True) 
        if res != None:
            t, phi, r = res 
            return self._trivial_factorization_with_n_phi(N, phi)

"""
Some composites
"""        
# B1, B2 = 2,0
Ns0 = [11861,557951,581851,1343807,
      3031993,4192681,5502053,6654671,
      12399797,14159471,16364501,20211713,
      22828523,44335457,50632823,57635863]

# B1,B2 = 10**5, 2
Ns1 = [384237701921,901500973759,6673435981363,
      11882099612003,15916520600381,17536455849431,
      18960823962761,21554451067267,33241863073423,
      55076499371497,57744632372831,67165261388497] 
      
Ns2 = [68072569242511,69684586201261,87756536366813,
      107458854058613,140967368218669,144086313393859,
      148077481187381,159872954386241,167863367550997,
      173112365219141,199390622342239,235255553137067,
      240522164849797,255119369061203,260436416434589,
      284911570131079,284984450489831,285341598723821,
      317418445093391,317554392679033,323219479761449,
      343840350739729,375275396864183,411429562199009,
      459621830953379,525220163614031]  
'''
RSA-59 RSA-79 RSA-99 RSA-119 RSA-100
'''
RSA = [71641520761751435455133616475667090434063332228247871795429,
       7293469445285646172092483905177589838606665884410340391954917800303813280275279,
       256724393281137036243618548169692747168133997830674574560564321074494892576105743931776484232708881,
       55519750778717908277109380212290093527311630068956900635648324635249028602517209502369255464843035183207399415841257091,
       1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139]         
#p1_list = [786766447,16375977287,81465486209,90824225567,862404698273,10452041068841,12697022389549,87996844075109,102010159808071,654472677526847,714033326215093,13051559264369500,13152735237439093,85817923293837151,131912444345458000]
#p2_list = [673815403,130260073,10937527,15171889,988483,109471,113489,11863,16141,919,631,83,67,13,11,]

def test(Ns,B2=0):
    Fib = Fibonacci()
    #times = []
    l = len(Ns)
    ff = 0
    n = 1
    tti = time.time()
    for N in Ns:
        if int(log2(N)) < 38:
            continue
        ti = time.time()
        B1, B2 = pow(10, int((log10(N)) // 2)-0), B2
        print("Test: %d, N: %d, log2(N): %d, B1: %d, B2: %d" % (n, N,int(log2(N)),B1,B2))
        P = Fib.factorization(N,B1,B2)
        if P != None:
            phi = (P[0]-1) * (P[1]-1)
            print("phi(N): %d" % phi)
            print("factors: %s" % str(P))
            ff += 1
        td = time.time() - ti
        ttd = time.time() - tti
        print("Runtime: %f\nFully factored:%d of %d" % (td,ff,l))
        print("Total Runtime: %f" % (ttd))
        n += 1
        print('------------------------------------------------------------------')


def test2():
  Fib = Fibonacci()
  N = 384237701921
  B1s = [10**x for x in range(6,3,-1)]
  B2s = [0,2,4,6,8]
  n=1
  tti = time.time()
  for B1 in B1s:
      for B2 in B2s:
          ti = time.time()
          print("Test: %d, N: %d, log2(N): %d, B1: %d, B2: %d" % (n, N,int(log2(N)),B1,B2))
          P = Fib.factorization(N,B1,B2)
          if P != None:
              phi = (P[0]-1) * (P[1]-1)
              print("phi(N): %d" % phi)
              print("factors: %s" % str(P))
              ff += 1
          td = time.time() - ti
          ttd = time.time() - tti
          print("Runtime: %f\nFully factored:%d of %d" % (td,ff,l))
          print("Runtime total: %f" % (ttd))
          n += 1


def test3(N, B2 = 0):
  Fib = Fibonacci()
  ti = time.time()
  B1, B2 = pow(10, int((log10(N)) // 2)), B2
  print("Test: N: %d, log2(N): %d, B1: %d, B2: %d" % (N,int(log2(N)),B1,B2))
  P = Fib.factorization(N,B1,B2)
  print(P)
  td = time.time() - ti
  print("Runtime: %f" % (td))

def test4(l,B2=0):
    Fib = Fibonacci(multitasked = True)
    #times = []
    ff = 0
    n = 1
    tti = time.time()
    for n in range(10,l):
        ti = time.time()
        n2 = 2**(n//2)
        n21 = (2**(n//2)-1)
        a = random.randint(n21,n2)
        p = next_prime(a)
        q = p
        while q == p:
           b = random.randint(n21,n2)
           q = next_prime(next_prime(b))
        N = p * q
        B1, B2 = pow(10, int((log10(N)) // 2)-0), B2
        print("Test: %d, N: %d, log2(N): %d, B1: %d, B2: %d" % (n, N,int(log2(N)),B1,B2))
        P = Fib.factorization(N,B1,B2)
        if P != None:
            phi = (P[0]-1) * (P[1]-1)
            print("phi(N): %d" % phi)
            print("factors: %s" % str(P))
            ff += 1
        td = time.time() - ti
        ttd = time.time() - tti
        print("Runtime: %f\nFully factored:%d of %d" % (td,ff,l))
        print("Total Runtime: %f" % (ttd))
        n += 1
        print('------------------------------------------------------------------')



if __name__=='__main__':
   #test(Ns0+Ns1+Ns2,B2=int(sys.argv[1]))
   #test2()
   #test3(int(sys.argv[1]))
   test4(100)
