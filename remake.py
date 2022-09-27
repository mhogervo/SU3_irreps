import sympy as sp

def poch(x,n):
    '''Return the Pochhammer symbol (x)_n, implemented as a recursion.'''
    if isinstance(n,int) and n >= 0:
        if n == 0: return 1
        elif n > 0: return (x+n-1)*poch(x,n-1)
    else:
        raise ValueError("n = {} must be a positive integer.".format(n))

def hermiteInnerProduct(s, t):
    '''
    If s = [i_1, i_2, ...] and t = [j_1, j_2, ...]  are two tuples that encode monomials in n coordinates,
    (such that s = x_1^{i_1} ... x_n^{i_n} etc.), return the L^2 inner product
    
    < s,t > = 1/\pi^{n/2} \int e^{-x^2} s(x) t(x)

    (the Hermite polynomials are orthonormal w.r.t. this inner product).
    '''
    if len(s) == len(t):
        sl, tl = list(s), list(t)
        out = 1
        while len(sl) > 0:
            n = sl.pop() + tl.pop()
            if n % 2 == 0:
                out *= poch(sp.Rational(1/2), n//2)
            else:
                return 0
        return out
    else:
        raise ValueError("Tuples {} and {} have different lengths.".format(s,t))

def buildGramMatrix(basis, innerProductFunction):
    '''Given a basis an an inner product function, return the Gram matrix.'''
    return [[innerProductFunction(x,y) for x in basis] for y in basis]

def addList(l):
    '''Take a tuple (k,l,m,n,...) and return a set containing (k+1,l,...), (k,l+1,m,...) etc.'''
    out = set()
    for i in range(len(l)):
        ln = list(l)
        ln[i] += 1
        out.add(tuple(ln))
    return out

def genBasis(ord,n=8):
    '''Return a set of all integer tuples of length n (i_1, ..., i_n) such that the sum of the i_k sum up to ord.'''
    if ord == 0:
        state = tuple([0 for _ in range(n)])
        return set([state])
    elif isinstance(ord,int) and ord >= 1:
        out = set()
        oldPartitions = genBasis(ord-1,n)
        for state in oldPartitions:
            out.update(addList(state))
        return out
    else:
        # if input is not sane 
        pass
    
def gramOfBasis(basis, gram):
    '''
    Let basis be a matrix of size (i,n) and gram a Gram matrix of size (n,n).
    Return the (i,i) matrix of inner products < basis_i, basis_j >.
    '''
    return basis.conjugate() * gram * basis.transpose()

def renormalizeMatrix(basis, gram):
    '''
    Given a matrix B = basis of size (i,n) and a gram matrix of size (n,n),
    return a new basis M' of size (i,n) s.t. (M')^* . gram . M^T = 1.
    '''
    if basis.cols == gram.rows == gram.cols:
        chol = gramOfBasis(basis, gram).cholesky()
        return chol.inverse() * basis
    else:
        raise ValueError("Size mismatch.")

def polyCoef(i,j,m,n):
    '''Coefficient of x^m y^n in the expansion (x+iy)^i (x-iy)^j.''' 
    for k in [i,j,m,n]:
        if (not isinstance(k,int)) or k < 0:
            raise ValueError("{} needs to be an integer >= 0.".format(k))
    if i+j == m+n:
        pref = pow(-1,j)*pow(sp.I,i+j-m)
        rest = sum([pow(-1,r)*sp.binomial(i,m-r)*sp.binomial(j,r) for r in range(m-i,j+1)])
        return pref*rest
    else:
        return 0
        
def expandPolyEven(s, t):
    '''
    Here s = [i_1, ...] encodes a "good monomial" of the form (x+iy)^{i_1} (x-iy)^{i_2}...
    and t = [n_1, ...] encodes a "basic monomial" x^{n_1} y^{n_2}..., 
    and it computes the coefficient in front of t when we expand s in terms of basic monomials.
    '''
    if len(s) % 2 == 1:
        raise ValueError("The << good monomial >> should have an even number of powers.")
    elif len(s) != len(t):
        raise ValueError("Both sets of polynomials should have the same number of variables.")
    else:
        if sum(s) == sum(t):
            sl, tl = list(s), list(t)
            out = 1
            while len(sl) > 0 and out != 0:
                i,j = sl.pop(0), sl.pop(0)
                m,n = tl.pop(0), tl.pop(0)
                out *= polyCoef(i,j,m,n)
            return out
        else:
            # if sum(s) != sum(t), the polynomials have different degree
            return 0

def expandPolyOdd(s, t):
    '''
    Here s = [i_1, ...] encodes a "good monomial" of the form x^{i_1} (y+iz)^{i_2} (y-iz)^{i_3}...
    and t = [n_1, ...] encodes a "basic monomial" x^{n_1} y^{n_2} z^{n_3}..., 
    and it computes the coefficient in front of t when we expand s in terms of basic monomials.
    '''
    if len(s) % 2 == 0:
        raise ValueError("The << good monomial >> should have an odd number of powers.")
    elif len(s) != len(t):
        raise ValueError("Both sets of polynomials should have the same number of variables.")
    else:
        if sum(s) == sum(t):
            sl, tl = list(s), list(t)
            a,b = sl.pop(0), tl.pop(0)
            if a==b: #the leading exponents need to be identical, otherwise it's 0
                out = 1
                while len(sl) > 0 and out != 0:
                    i,j = sl.pop(0), sl.pop(0)
                    m,n = tl.pop(0), tl.pop(0)
                    out *= polyCoef(i,j,m,n)
                return out
            else:
                return 0
        else:
            # if sum(s) != sum(t), the polynomials have different degree
            return 0

#class SU3poly
# lev = 3

# bas = tuple(genBasis(lev))
# gram = sp.ImmutableSparseMatrix(buildGramMatrix(bas, innerProduct))
# dimBas = len(bas)
# print(dimBas)
# # construct a new basis of states that have definite isospin and hypercharge quantum numbers
# definiteBasis = sp.ImmutableSparseMatrix([[expandPoly(s,t) for s in bas] for t in bas]).transpose()
# # orthonormalize w.r.t. the gram matrix
# definiteBasis = renormalizeMatrix(definiteBasis, gram)

# # check that the new basis is properly normalized
# gramOfBasis(definiteBasis, gram) == sp.eye(dimBas)

