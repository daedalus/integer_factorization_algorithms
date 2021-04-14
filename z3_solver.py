from z3 import Solver, Int ,set_param
from gmpy2 import isqrt
set_param('parallel.enable', True)

def z3_solve(self, n, timeout_amount):
        s = Solver()
        s.set("timeout", timeout_amount * 1000)
        p = Int("x")
        q = Int("y")
        i = int(isqrt(n))
        if i**2 == n: # check if we are dealing with a perfect square otherwise try to SMT.
            return i,i
        s.add(p * q == n, p > 1, q > i, q > p) # In every composite n=pq,there exists a p>sqrt(n) and q<sqrt(n).
        try:
            s_check_output = s.check()
            res = s.model()
            return res[p].as_long(), res[q].as_long()
        except:
            return None, None
            
print(z3_solve(int(sys.argv[1])))
