def Analytical_Solution(u_exact, v_exact, p_exact, xu, yu, xv, yv, xp, yp, T):
    '''Calculates the analytical solution of Taylor Green Vortex at time T.
    
    Params:
    -------
    u_exact, v_exact, p_exact   2D array of float, exact solution grids
    xu, yu                      2D array of float, u velocity grid points
    xv, yv                      2D array of float, v velocity grid points
    xp, yp                      2D array of float, pressure grid points
    T                           float, solution Time
    
    Returns:
    --------
    u, v                        2D array of float, analytical u and v at T
    p                           2D array of float, analytical p at T
    '''
    u = np.empty_like(u_exact, dtype=float)
    v = np.empty_like(v_exact, dtype=float)
    p = np.empty_like(p_exact, dtype=float)


    for i in range(0,Ny+2):
        for j in range(0,Nx+3):
            u[i,j] = -np.exp(-2.0*T)*np.cos(xu[j])*np.sin(yu[i])



    for i in range(0,Ny+3):
        for j in range(0,Nx+2):
            v[i,j] = np.exp(-2.0*T)*np.sin(xv[j])*np.cos(yv[i])


    for i in range(0,Ny):
        for j in range(0,Nx):
            p[i,j] = -0.25*np.exp(-4.0*T)*(np.cos(2.0*xp[j]) +\
                np.cos(2.0*yp[i]))
            
    return u, v, p