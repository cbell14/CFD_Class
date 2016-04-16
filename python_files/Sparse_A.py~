def Build_Sparse_A_Neumann(Nx, Ny):
    """Builds the A matrix for a 5 point stencil Poisson Solver
    
    Params:
    -------
    Nx, Ny     float, number of points in x and y
    
    Returns:
    -------
    A          Penta-diagonal sparse matrix A
    """
    
    A = sparse.diags([-4]*Nx*Ny, 0) #set leading Diaganol
    A = A + sparse.diags([1]*(Nx-1)+([0]+[1]*(Nx-1))*(Ny-1), 1) #set first diagonal above main
    A = A + sparse.diags([1]*(Nx-1)+([0]+[1]*(Nx-1))*(Ny-1), -1) #set first diagonal below main
    A = A + sparse.diags([1]*(Nx*Ny-Nx), Nx) #sets second diagonal above main
    A = A + sparse.diags([1]*(Nx*Ny-Nx), -Nx) #sets second diagonal below main
    
    return A