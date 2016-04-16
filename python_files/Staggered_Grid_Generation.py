def Generate_Grid(Lx, Ly, Nx, Ny, x_s, y_s):
    """Generates a 2D staggered Grid
    
    Params:
    -------
    Lx, Ly     float, length of grid in x and y
    Nx, Ny     float, number of points in x and y not includeing ghost cells
    xs, ys     float, starting points in x and y
    
    Returns:
    -------
    dx, dy     float, step size in x and y
    x, y       1D array of float, x and y points
    xp, yp     1D array of float, pressure point grid in x and y
    p_exact    2D array of float, array for exact pressure
    xu, yu     1D array of float, u velocity grid in x and y
    u_exact    2D array of float, array for exact u
    ut         2D array of float, array for time varying u
    xv, yv     1D array of float, v velocity grid in x and y
    v_exact    2D array of float, array for exact v
    vt         2D array of float, array for time varying v
    """
    
    dx = Lx/(Nx)
    dy = Ly/(Ny)
    
    x = np.zeros((Nx+2,1), dtype=float)
    y = np.zeros((Ny+2,1), dtype=float)
        
    # Pressure Points
    xp = np.zeros(Nx, dtype=float)
    yp = np.zeros(Ny, dtype=float)
    p_exact = np.zeros((Ny,Nx), dtype=float)


    # u velocity points
    xu = np.zeros(Nx+3, dtype=float)
    yu = np.zeros(Ny+2, dtype=float)
    u_exact = np.zeros((Ny+2,Nx+3), dtype=float)
    #ut = np.zeros_like((Ny+2,Nx+3,nt), dtype=float)

    # v velocity points
    xv = np.zeros(Nx+2, dtype=float)
    yv = np.zeros(Ny+3, dtype=float)
    v_exact = np.zeros((Ny+3,Nx+2), dtype=float)
    #vt = np.zeros_like((Ny+3,Nx+2,nt), dtype=float)
    
    #Pressure
    #for i in range(0,Nx+2):
    xp = np.linspace(dx/2.0,Lx-dx/2.0,Nx)
    yp = np.linspace(dy/2.0,Ly-dy/2.0,Ny)

    #u Velocity
    xu = np.linspace(-dx,Lx+dx,Nx+3)
    yu = np.linspace(-dy/2.0,Ly+dy/2.0,Ny+2)

    #v Velocity
    xv = np.linspace(-dx/2.0,Lx+dx/2.0,Nx+2)
    yv = np.linspace(-dy,Ly+dy,Ny+3)
    
    return x, y, xp, yp, p_exact, xu, yu, u_exact, xv, yv, v_exact, dx, dy