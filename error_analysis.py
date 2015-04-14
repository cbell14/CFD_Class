import nump as np

def error_analysis(Nx,Ny,x_vals,exact_sol):
    """Computes L2 norm of calculated solution
    
    Params:
    ------
    Nx, Ny      Number of points
    x_vals      calculated values
    exact_sol   Analytical Solution
    
    Returns:
    -------
    num         Number of total points
    normdiff    L2 norm of calculated solution
    """

    num = Nx * Ny
    temp = np.reshape(x_vals, (Ny,Nx))
    diff = exact_sol - temp
    normdiff = np.linalg.norm(diff,ord=2)
    
    return num, normdiff