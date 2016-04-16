def Poisson_Solver_Neumann(Nx, Ny):
    """Solves the 2D Poisson equation implicitly on a staggered grid
    using Neumann Boundary Conditions
    
    Params:
    ------
    u, v           2D array of float, x and y velocities
    Nx, Ny         float, Number of segments in x and y
    dt             float, time step size
    T              float, current time
    X, Y           2D array of float, meshgrid
    
    Returns:
    -------
    ANeum        2D array of float, A matrix with Neumann conditions
    f_RHSn       1D array of float, f(x,y) for Neumann conditions
    """
    #Building A
    ANeum = np.zeros((Nx*Ny,Nx*Ny),dtype=float)

    a = 1.0
    b = 1.0
    c_int = -4.0 
    c_edge = -3.0
    c_corner = -2.0
    d = 1.0
    e = 1.0

    #Set corner points
    ANeum[0,0] = c_edge
    ANeum[-1,-1] = c_corner
    ANeum[Nx*Ny-Ny,Nx*Ny-Ny] = c_corner
    ANeum[Ny-1,Ny-1] = c_corner

    #Set edges in first block
    for j in range(1,Ny-1):
        ANeum[j,j] = c_edge
        j +=j

    #Set edges in last block
    for j in range((Nx*Ny)-Ny,(Nx*Ny)-2):
        ANeum[j+1,j+1] = c_edge
        j +=j

    #Set edges along main diagonal except for first block
    for j in range(Nx+1,Ny*Nx):
        if j % Nx ==0:
            ANeum[j-1,j-1] = c_edge
        j +=j

    #Set edges on main diagonal except for last block
    for j in range(Nx,(Ny*Nx)-Nx):
        if j % Nx ==0:
            ANeum[j,j] = c_edge
        j +=j

    #Second diagonal above and below diagonal
    for j in range(Ny,Ny*Nx):
        ANeum[j,j-Ny] = a
        ANeum[j-Nx,j] = e
        j +=j

    #first diagonal below main diagonal   
    for j in range(1,Ny*Nx):
        if j % Ny ==0:
            ANeum[j,j-1] = 0
        else:
            ANeum[j,j-1] = b
        j +=j

    #first diagonal above main diagonal
    for j in range(0,Ny*Nx):
        if j % Nx ==0:
            ANeum[j-1,j] = 0
        else:
            ANeum[j-1,j] = d
        j +=j

    #Main Diagonal
    for j in range(0,Ny*Nx):
        if ANeum[j,j] ==0:
            ANeum[j,j] = c_int
    
    return ANeum