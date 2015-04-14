import numpy as np

def Conj_Grad(A,b,max_iter):
    """Solves Ax=b using the conjugate gradient method
    
    Params:
    ------
    A         Positive Definite and symmetric coefficient matrix m by m
    b         RHS of matrix m by 1
    max_iter  maximum number of iterations before breaking
    
    Returns:
    ------
    x   
    """
    x = np.zeros((A.shape[0],1), dtype=float)
    r = b.copy()
    rtr = np.dot(r.T,r)
    d = r.copy()
    k = 0
    error = 100
    
    while abs(error) > 1e-4:
        if  k > max_iter:
            print("Maximum number of iterations reached for Conjugate Gradient")
            break
        else:
            Ad = np.dot(A,d)
            alpha = rtr / np.dot(d.T,Ad)
            x = x + alpha * d
            r = r - alpha * Ad
            rtrold = rtr
            rtr = np.dot(r.T,r)
            beta = rtr / rtrold
            d = r + beta * d

            error = np.sqrt(rtr)/np.linalg.norm(b, ord=2)
            k +=k

    bnorm = np.linalg.norm(b,ord=2)
    total_error = np.linalg.norm(r,ord=2) / bnorm
    return x, total_error