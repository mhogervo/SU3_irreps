def checkPosInt(n,n0=0):
    if (not isinstance(n,int)):
        raise TypeError(f"{n} is not an integer.")
    elif n < n0:
        raise ValueError(f"{n} is an integer < {n0}")
    else:
        pass

def shift(ar,k=1):
    '''
    Take a list/tuple and shift all entries to the right k times.
    k can be any integer (including negative).
    '''
    if not isinstance(k,int):
        raise TypeError(f"Input {k} is not an integer.")
    elif not (isinstance(ar,list) or isinstance(ar,tuple)):
        raise TypeError(f"Input {ar} should be a list/tuple.")
    else:
        pass

    a, L = list(ar), len(ar)
    if k==0 or L <= 1:
        return ar
    elif k==1:
        n = a.pop()
        a = [n] + a
        if isinstance(ar,tuple):
            a = tuple(a)
        return a
    elif 1 < k < L:
        return shift(shift(ar,1),k-1)
    else:
        return shift(ar, k % L)        
    
def flip(ar):
    '''
    Take a list/tuple and reverse it.
    '''
    if isinstance(ar,tuple):
        return tuple(reversed(ar))
    elif isinstance(ar,list):
        return list(reversed(ar))
    else:
        raise TypeError(f"{ar} is not a tuple/list.")

def orbit(ar,bracelet=False):
    '''
    Given an array, generate the orbit.
    '''
    L = len(ar)
    cyclic_orbit = {tuple(shift(ar,k)) for k in range(L)}
    if bracelet:
        flipped = {flip(x) for x in cyclic_orbit}
        cyclic_orbit.update(flipped)
    return frozenset(cyclic_orbit)

def orbitRep(ar,bracelet=False):
    '''
    Given a array, generate the orbit, and pick a representative.
    '''
    orbitOfAr = sorted(list(orbit(ar,bracelet)))
    return orbitOfAr[0]

def allNecklaces(n,k=2,bracelet=False):
    '''
    Generate all necklaces with k beads of length n.
    '''
    checkPosInt(n), checkPosInt(k,1)

    if n==0:
        return [()] #empty necklace
    elif n==1:
        return [(i,) for i in range(k)]
    elif n>1:
        out = set()
        shortNecklaces = allNecklaces(n-1,k,bracelet)
        for x in shortNecklaces:
            orb = orbit(x,bracelet)
            new = {y + (i,) for y in orb for i in range(k)}
            for y in new:
                out.add(orbitRep(y,bracelet))
        return sorted(list(out))
    
def allBracelets(n,k=2):
    return allNecklaces(n,k,bracelet=True)

###############################################################
###############################################################
####
#### For SU(N) matrices, there is a relation
####
####     X^N = 1/N Tr(X^N) + lower traces
####
#### which means that many necklaces don't contribute. We need a functions
#### that only keeps appropriate necklaces.
####
###############################################################
###############################################################

def largestDomain(necklace):
    # A necklace is of the form (i1, i2, i3, ...). A domain is a group of neighboring beads of the same color.
    # This function computes the size of the largest domain.
    
    # the array walls contains the positions of all the domain walls
    walls, L = [], len(necklace)
    for pos in range(L):
        i,j = necklace[pos-1], necklace[pos]
        if i != j: walls.append(pos)
    
    if walls == []:
        return len(necklace)
    else:
        number_walls = len(walls)
        sizeOfDomains = [((walls[i] - walls[i-1]) % L) for i in range(number_walls)]
        return max(sizeOfDomains)    

def validSingleTraces(n,k,N):
    '''
    Generate all necklaces with k beads of length n, and at most N-1 beads of the same color appear in a row.
    For n=N we make an exception, since (i,i,...,i) [n copies of letter i] is a valid single-trace operator.
    '''
    checkPosInt(n), checkPosInt(k,1), checkPosInt(N,2)
     
    if n==0:
        return [()] #empty necklace
    elif n==1:
        return [] #throw away tr(j) = 0
    elif n>1:
        # function that checks if a necklace does not contain
        # a redundant operator:
        check = lambda necklace : largestDomain(necklace) < N
        
        rawNecklaces = allNecklaces(n,k)
        out = set(filter(check, rawNecklaces))

        if n==N:
            constantNecklace = lambda j : tuple([j for _ in range(N)])
            out.update(set([constantNecklace(j) for j in range(k)]))

        return sorted(list(out))
    
###############################################################
###############################################################
####
#### This last part contains number-theoretical functions that can be used
#### to check the results.
####
###############################################################
###############################################################

def gcd(a,b):
    # from Knuth, Euclidean algorithm
    checkPosInt(a,1), checkPosInt(b,1)
    while b != 0:
        a, b = b, a % b
    return a

def eulerPhi(n):
    '''
    The Euler totient function, its primary definition.
    '''
    coprime = [m for m in range(1,n+1) if gcd(m,n) == 1]
    return len(coprime)

def divisors(n):
    '''
    Set of all divisors of a positive integer n.
    '''
    checkPosInt(n,1)
    return [m for m in range(1,n+1) if n % m == 0]

def necklacePoly(n,k=2):
    '''
    Number-theoretical counting for the number of necklaces of length n
    made with k beads. See Wikipedia.
    '''
    checkPosInt(n), checkPosInt(k,1)
    if n == 0:
        return 1
    else: # n>=1
        return sum([eulerPhi(d)*pow(k,n//d) for d in divisors(n)])//n
        

def braceletPoly(n,k=2):
    '''
    Number-theoretical counting for the number of bracelets of length n
    made with k beads. See Wikipedia.
    '''
    if n == 0:
        return 1
    elif n % 2 == 0:
        return (necklacePoly(n,k) + ((k+1)*pow(k,n//2))//2)//2
    else: # n % 2 == 1:
        return (necklacePoly(n,k) + pow(k,(n+1)//2))//2


