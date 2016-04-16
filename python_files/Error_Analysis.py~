def error_analysis(num_sol,exact_sol):
    """Computes L2 norm of calculated solution
    
    Params:
    ------
    N_cx, N_cy  Number of cell centered points
    x_vals      calculated x_values
    exact_sol   Analytical Solution
    
    Returns:
    -------
    num         Number of total points
    normdiff    L2 norm of calculated solution
    """
    diff = num_sol - exact_sol
    norm2diff = np.linalg.norm(diff,ord=2)/np.sqrt(np.sum((exact_sol**2)))
    normdiff = np.linalg.norm(diff,ord=1)/np.sqrt(np.sum((exact_sol**2)))
    
    return normdiff, norm2diff