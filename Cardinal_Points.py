def Cardinal_Points(x, Nx, Ny):
    """Takes input and finds all north, south, east,  west, and central points
    
    Params:
    ---------
    x             2D array of float, target vector
    Nx, Ny   float, number of points in x and y
    
    Returns:
    ---------
    x_n      2D array of float
    x_s      2D array of float
    x_e      2D array of float
    x_w     2D array of float
    x_p      2D array of float
    """
    
    x_n = np.zeros((Nx, Ny), dtype=float)
    x_s = np.zeros((Nx, Ny), dtype=float)
    x_e = np.zeros((Nx, Ny), dtype=float)
    x_w = np.zeros((Nx, Ny), dtype=float)
    x_p = np.zeros((Nx, Ny), dtype=float)
    
    for i in range(0,Nx-1):
        for j in range(0,Ny):
            x_n[i,j] = x[i+1,j]
            
    for i in range(1,Nx):
        for j in range(0,Ny):
            x_s[i,j] = x[i-1,j]
            
    for i in range(0,Nx):
        for j in range(0,Ny-1):
            x_e[i,j] = x[i,j+1]
            
    for i in range(0,Nx):
        for j in range(1,Ny):
            x_w[i,j] = x[i,j-1]
    for i in range(0,Nx):
        for j in range(0,Ny):
            x_p[i,j] = x[i,j]
        
    return x_n, x_s, x_e, x_w, x_p