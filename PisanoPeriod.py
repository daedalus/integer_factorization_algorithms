#!/usr/bin/env python3
"""
Integer factorization with pisano period
Original repo https://github.com/wuliangshun/IntegerFactorizationWithPisanoPeriod/
White paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8901977
"""

import random
import time
from gmpy2 import isqrt, gcd
import sys
sys.setrecursionlimit(5000)

class Fibonacci:
    def __init__(self):
        pass
    

    def _fib_res(self,n,p):
    	if n == 0:
            return (0, 1)
    	else:
            a, b = self._fib_res(n >> 1,p)
            c = ((a % p) * ((b << 1) - a) % p) % p
            d = (pow(a, 2, p) + pow(b, 2, p)) % p
            if n & 1 == 0: 
                return (c, d)
            else:
                return (d, c + d)
    

    def get_n_mod_d(self,n,d):
        if n < 0:
            ValueError("Negative arguments not implemented")
        return self._fib_res(n,d)[0]
    
    
    def binary_search(self,L,n):
        left = 0
        right = len(L) - 1
        while left <= right:
            mid=(left + right) >> 1
            if n == L[mid]:
                return mid
            elif n < L[mid]:
                right = mid - 1
            else:
                left = mid + 1
        return -1
    
    
    def sort_list(self,L):
        from operator import itemgetter
        indices, L_sorted = zip(*sorted(enumerate(L), key=itemgetter(1)))
        return list(L_sorted),list(indices)
   
 
    def get_period_bigint(self, N, min_accept, xdiff, verbose = False):            
        search_len = int(pow(N, (1.0 / 6) / 100))
        
        if search_len < min_accept:
            search_len = min_accept
  
        if verbose:
            print('search_len:%d'%(search_len))
        
        starttime = time.time()
        diff = xdiff 
        p_len = int((len(str(N)) + diff) >> 1) + 1
        begin = N - int('9'*p_len) 
        if begin <= 0:
            begin = 1
        end = N + int('9' * p_len)
    
        if verbose:    
            print('search begin:%d,end:%d'%(begin, end))
                
        rs = [self.get_n_mod_d(x, N) for x in range(search_len)]
        rs_sort, rs_indices = self.sort_list(rs) 
             
        if verbose:    
            print('sort complete! time used:%f secs' % (time.time() - starttime))
                
        T = 0
        has_checked_list = []

        while True:       
            randi = random.randint(begin,end)            
            if (not randi in has_checked_list):
                has_checked_list.append(randi)
            
                res = self.get_n_mod_d(randi, N)
                inx = self.binary_search(rs_sort, res)
                                    
                if inx > -1:                
                    res_n = rs_indices[inx]
                    T = randi - res_n
                     
                    if self.get_n_mod_d(T, N) == 0:
                        t = time.time() - starttime
                      
                        if verbose:
                            print('For N = %d Found T:%d, randi: %d, time used %f secs.' % (N , T, t, randi))
                        return t, T, randi
                    else:
                        if verbose:
                            print('For N = %d\n Found res: %d, inx: %d, res_n: %d , T: %d\n but failed!' % (N, res, inx, res_n, T))
                        #return None

  
    def _trivial_factorization_with_n_t(self, N, T):
        M = N - T + 1
        d = N
        p1 = []

        M2 = pow(M,2)
        M2p4d = M2 + 4*d
        M2m4d = M2 - 4*d

        if M2m4d > 0:
            iM2m4d = isqrt(M2m4d)
            p1.append((M + iM2m4d) >> 1)
            p1.append((M - iM2m4d) >> 1)

        if M2p4d > 0:
            iM2p4d = isqrt(M2p4d)
            p1.append((M + iM2p4d) >> 1)
            p1.append((M - iM2p4d) >> 1)

        for p in p1:
            g = gcd(p,N)
            if N > g > 1:
                return g,N//g
   

    def factorization(self, N, min_accept, xdiff):
        res = self.get_period_bigint(N, min_accept, xdiff, verbose=True) 
        if res != None:
            t, T, r = res
            return self._trivial_factorization_with_n_t(N, T)

"""
Some composites
"""        
Ns = [11861,557951,581851,1343807,
      3031993,4192681,5502053,6654671,
      12399797,14159471,16364501,20211713,
      22828523,44335457,50632823,57635863,
      384237701921,901500973759,6673435981363,
      11882099612003,15916520600381,17536455849431,
      18960823962761,21554451067267,33241863073423,
      55076499371497,57744632372831,67165261388497,
      68072569242511,69684586201261,87756536366813,
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
    
if __name__=='__main__':
    fib = Fibonacci()
    times = []
    for N in (Ns + RSA):
        print("N: %d" % N)
        P = fib.factorization(N,2,0)
        if P != None:
            phi = (P[0]-1) * (P[1]-1)
            print("phi(N): %d" % phi)
            print("factors: %s" % str(P))


        print('---------------------------------')