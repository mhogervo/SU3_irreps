from sympy import ImmutableSparseMatrix, I, sqrt
import sympy as sp

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

def commutator(A,B):
    return A*B - B*A

def GellMann():
    '''For N=3, spit out a basis of generators T_a and structure constants f_{abc}.'''
    global N
    if N != 3: raise ValueError("Only works for SU(3), not SU({})".format(N))
   
    basis = []
    M  = ImmutableSparseMatrix(N,N, {(0,0): 1, (1,1): -1})
    basis.append(M)
    M  = ImmutableSparseMatrix(N,N, {(0,0): 1, (1,1): 1, (2,2): -2})
    basis.append(M)
    for i in range(N):
        for j in range(i):
            M = ImmutableSparseMatrix(N,N, {(i,j): 1, (j,i): 1})
            basis.append(M)
            M = ImmutableSparseMatrix(N,N, {(i,j): I, (j,i): -I})
            basis.append(M)
    basis = [M/M.norm() for M in basis]
    sc = [[[ (A*commutator(B,C)).trace() for A in basis] for B in basis] for C in basis ]
    sc = I*sp.ImmutableSparseNDimArray(sc)
    return basis, sc

def action(a,psi,lu):
    '''Action of the a'th generator on the state psi. The result is looked up
    using a lookup table lu.'''
    global dimAdjoint, structureConstants
    out = sp.zeros(1,len(lu))
    for b in range(dimAdjoint):
        for c in range(dimAdjoint):
            m,n = psi[b], psi[c]
            if n > 0 and b != c:
                psiNew = list(psi)
                psiNew[b] += 1
                psiNew[c] -= 1
                psiNew = tuple(psiNew)
                pref = I*structureConstants[a,b,c]*sp.sqrt((m+1)*n)
                out[lu[psiNew]] += pref
    return out

def isUnitary(M):
    return M * M.adjoint() == sp.eye(M.rows)

def genSUNmatrix(a,lu):
    '''The action of the a'th generator of SU(N) represented on a basis indicated by a lookup table lu.'''
    nn = len(lu)
    M = sp.MutableSparseMatrix(nn,nn,{})
    stateList = list(lu.keys())
    for i in range(nn):
        state = stateList[i]
        M[:,i] = [x for x in action(a,state,lu)]
    return sp.ImmutableSparseMatrix(M)

def normVectorsSymPy(vecs):
    '''Given a set of vectors in SymPy format, return an orthonormalized matrix.'''
    return sp.Matrix(sp.GramSchmidt([x.transpose() for x in vecs],True))

def generateRep(ord):
    # first make a list of all allowed polynomials, label them and store as a dict/lookup table
    polyBasis = list(genBasis(ord,dimAdjoint))
    polyBasis = sorted(polyBasis)
    lookup = {}
    for n, state in enumerate(polyBasis):
        lookup[state] = n
    nn = len(lookup)
    
    hermiteMatrices = [genSUNmatrix(a,lookup) for a in range(dimAdjoint)]
    
    isospinMatrix = sp.sqrt(2) * hermiteMatrices[0]
    hyperchargeMatrix = sp.sqrt(6)/3 * hermiteMatrices[1]
    casimirMatrix = sp.zeros(nn,nn)
    for M in hermiteMatrices: casimirMatrix += M*M
    
    ###
    ### compute the change-of-basis matrix that diagonalizes isospin,
    ### hypercharge and the Casimir
    ###
    
    ### fix the order of operations:
    matrixList = (casimirMatrix, hyperchargeMatrix, isospinMatrix)
    return isospinMatrix
    U = sp.ImmutableSparseMatrix(0,nn,{})
    print("Diagonalizing...")
    initSpectrum = matrixList[0].eigenvects()       ## this is computationally
                                                    ## the only time-consuming job
    print("...done.")
    for x in initSpectrum:
        _, _, spec1 = x
        spec1 = normVectorsSymPy(spec1)
        block1 = spec1 * matrixList[1] * spec1.adjoint()
        for xx in block1.eigenvects():
            _, _, spec2 = xx
            spec2 = normVectorsSymPy(spec2)*spec1
            #print("iso: ",isUnitary(specIso),specIso.shape)
            block2 = spec2 * matrixList[2] * spec2.adjoint()
            for xxx in block2.eigenvects():
                _, _, spec3 = xxx
                spec3 = normVectorsSymPy(spec3)*spec2
                U = U.row_insert(-1,spec3)
      
    return lookup, hermiteMatrices, casimirMatrix, U

N = 3
dimAdjoint = pow(N,2) - 1
_, structureConstants = GellMann()
#print(basis)
#print(structureConstants)

ord=3
# lookup, hermiteRep, casimirMatrix, U = generateRep(ord)
M = generateRep(ord)

#nn = len(lookup)
#isospinMatrix = sp.sqrt(2) * hermiteRep[0]
#hyperchargeMatrix = sp.sqrt(6)/3 * hermiteRep[1]
#casimirMatrix = sp.zeros(nn,nn)
#for M in hermiteRep: casimirMatrix += M*M
#print(isUnitary(U),U.shape)
#print("There are {} polynomials in this basis.".format(len(lookup)))
# check that the new basis diagonalizes the matrices in question
#[ (U*M*U.adjoint()).is_diagonal() for M in (casimirMatrix, hermiteRep[0], hermiteRep[1]) ]

# check that the R_a obey [ R_a, R_b ] = - i f_{abc} R_c

def checkHermite(a,b):
    global N, dimAdjoint, hermiteRep
    nn = hermiteRep[0].rows
    out = sp.zeros(nn)
    for c in range(dimAdjoint):
        out += -I*structureConstants[a,b,c]*hermiteRep[c]
    return out

[[ (commutator(hermiteRep[a],hermiteRep[b]) - checkHermite(a,b)).norm() for a in range(dimAdjoint)] for b in range(dimAdjoint)]

newHermiteRep = [ U * M * U.adjoint() for M in hermiteRep ]
isn = U * isospinMatrix * U.adjoint()
hcn = U * hyperchargeMatrix * U.adjoint()
cas = U * casimirMatrix * U.adjoint()
print("Finished diagonalizing isospin, hypercharge and the Casimir.")