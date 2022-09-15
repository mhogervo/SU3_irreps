import sympy as sp

def buildGramMatrix(basis, innerProductFunction):
    '''
    Given a basis {x} (as an iterator) and an inner product function f(x,y) of two arguments, 
    return the Gram matrix < x | f | y >.
    '''
    mat = [[innerProductFunction(x,y) for x in basis] for y in basis]
    return sp.ImmutableSparseMatrix(mat)


def gramOfBasis(basis, gram):
    '''
    Let basis be a SymPy matrix of size (i,n) and gram a (hermitian) Gram matrix of size (n,n).
    Return the (i,i) matrix of inner products < basis_i, basis_j >.
    '''
    if basis.cols == gram.rows == gram.cols:
        return basis.conjugate() * gram * basis.transpose()
    else:
        raise ValueError("Matrix size mismatch.")

def renormalizeMatrix(basis, gram):
    '''
    Given a matrix B = basis of size (i,n) and a gram matrix of size (n,n),
    return a new basis M' of size (i,n) s.t. (M')^* . gram . M^T = 1.
    '''
    if basis.cols == gram.rows == gram.cols:
        chol = gramOfBasis(basis, gram).cholesky()
        return chol.inverse() * basis
    else:
        raise ValueError("Matrix size mismatch.")

def isUnitary(M):
    return M * M.adjoint() == sp.eye(M.rows)

def normVectors(vecs):
    '''Given a set of vectors in SymPy format, return an orthonormalized matrix.'''
    return sp.Matrix(sp.GramSchmidt([x.transpose() for x in vecs],True))


def simul_diag(matArray, vecs):
    ''' matArray should be a list of matrices of size (n,n); vecs should be a matrix of size (i,n).'''
    if len(matArray) == 1:
        M = matArray[0]
        M_restr = vecs.conjugate() * M * vecs.transpose()
        spec = normVectors(M_restr.eigenvects())
        return spec
        
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

