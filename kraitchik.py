from gmpy2 import *

def kraitchik(n):
    x = isqrt(n)
    while True:
        k = 1
        s = x * x - k * n
        while s >= 0:
            if is_square(s):
                y = isqrt(s)
                z = x + y
                w = x - y
                if z % n != 0 and w % n != 0:
                    return gcd(z,n),gcd(w,n)
            s = x * x - k * n
            k += 1
        x+=1

print(kraitchik(2041))
print(kraitchik(65537*65539))

