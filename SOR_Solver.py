import numpy as np

def SOR(A,b,omega,max_iter,sigma):
    """Solves Ax=b using the successive over relaxation method, if omega=1
    this methods becomes Gauss-Seidel
    
    Params:
    ------
    A         square coefficient matrix n by m
    b         RHS of matrix n by 1
    omega     relaxation parameter
    max_iter  maximum number of iterations before breaking
    sigma     desired convergence level
    
    Returns:
    ------
    x       solved variable matrix
    """
    #generating first guess matrix
    x0 = np.zeros((A.shape[0],1), dtype=float) 
    x1 = np.empty((A.shape[0],1), dtype=float)
    I = np.eye(A.shape[0],dtype=float)
    L = np.tril(A,-1)
    U = np.triu(A,1)
    D = np.diag(A)
    temp = np.linalg.inv(L + np.diagflat(D))
    q = np.dot(temp,b)
    k = 0
    error = 100.0
    
    while abs(error) > sigma:
        if k > max_iter:
            print("Maximum number of iterations reached for SOR")
            break
        else:
            x1 = np.dot((omega*(I - np.dot(temp,A))) + ((1-omega)*I),x0) + omega*q
            diff = x1 - x0
            error = np.linalg.norm(diff, ord=2) 
            x0 = x1.copy()
            k += 1
            
    bnorm = np.linalg.norm(b,ord=2)
    total_error = np.linalg.norm(diff,ord=2) / bnorm
    return x1, total_error