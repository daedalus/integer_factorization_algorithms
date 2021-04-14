# Code borrowed and adapted from the wikipedia: https://en.wikipedia.org/wiki/Shanks%27s_square_forms_factorization
# It may contain bugs

def SQUFOF(N):
    s = int( isqrt(N)+0.5)
    if (s**2 == N):
        return s
    for k in range(0,len(multiplier)):
        D = multiplier[k]*N
        Po = Pprev = P = isqrt(D)
        Qprev = 1
        Q = D - Po**2
        L = int(2 * isqrt(2*s))
        B = 3 * L
        for i in range(2,B):
            b = int((Po + P)//Q)
            P = b*Q - P
            q = Q
            Q = Qprev + b*(Pprev - P)
            r = int(isqrt(Q)+0.5)
            if (not(i & 1) and ((r**2) == Q)):
                break
            Qprev = q
            Pprev = P
        if (i >= B):
            continue
        b = ((Po - P)//r)
        Pprev = P = b*r + P
        Qprev = r
        Q = (D - Pprev**2)//Qprev
        i = 0
        c = True
        while(c):
            b = int((Po + P)//Q)
            Pprev = P
            P = b*Q - P
            q = Q;
            Q = Qprev + b*(Pprev - P)
            Qprev = q
            i+=1
            c = (P != Pprev)
        r = gcd(N, Qprev)
        if (1 < r < N):
            return r
    return -1

N=int(sys.argv[1])
print(SQUFOF(N))
@daedalus
 
