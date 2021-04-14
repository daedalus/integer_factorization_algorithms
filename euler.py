# /usr/bin/env python
# code taken from https://maths.dk/teaching/courses/math357-spring2016/projects/factorization.pdf

from gmpy2 import gcd, isqrt

def euler(self, n):
        if n % 2 == 0:
            return (n / 2, 2) if n > 2 else (2, 1)
        end = isqrt(n)
        a = 0
        solutionsFound = []
        firstb = -1
        while a < end and len(solutionsFound) < 2:
            bsquare = n - a ** 2
            if bsquare > 0:
                b = isqrt(bsquare)
                if (b ** 2 == bsquare) and (a != firstb) and (b != firstb):
                    firstb = b
                    solutionsFound.append([int(b), a])
            a += 1
        if len(solutionsFound) < 2:
            return -1
        a = solutionsFound[0][0]
        b = solutionsFound[0][1]
        c = solutionsFound[1][0]
        d = solutionsFound[1][1]
        k = gcd(a - c, d - b)
        h = gcd(a + c, d + b)
        m = gcd(a + c, d - b)
        l = gcd(a - c, d + b)
        n = (k ** 2 + h ** 2) * (l ** 2 + m ** 2)
        return [int(k ** 2 + h ** 2) // 2, int(l ** 2 + m ** 2) // 2]
        
print(euler(int(sys.argv[1])))
