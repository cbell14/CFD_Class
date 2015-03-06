import numpy as np

def Jacobi_Solver(A,b,max_iter,sigma):
    """Uses the jacobi method to solve Ax=b
    
    Params:
    ------
    A         square coefficient matrix n by m
    b         RHS of matrix n by 1
    max_iter  maximum number of iterations before breaking
    sigma     desired convergence level
    
    Returns:
    ------
    x       solved variable matrix
    """
    
    #create intial guess x the size of A's first column
    x0 = np.zeros((A.shape[0],1), dtype=float) 
    x1 = np.empty((A.shape[0],1), dtype=float) #ditto
    I = np.eye(A.shape[0])
    D = np.diag(A)
    q = np.dot(np.diagflat(1/D),b)
    k = 0
    error = 100.0
    
    while abs(error) > sigma:
        if k > max_iter:
            print("Maximum number of iterations reached")
            break
        else:
            x1 = np.dot((I - np.dot(np.diagflat(1/D),A)),x0) + q
            diff = x1 - x0
            error = np.linalg.norm(diff, ord=2) 
            x0 = x1.copy()
            k += 1
    return x1