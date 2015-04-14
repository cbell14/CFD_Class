def Vorticity(u, v, dx, dy):
    """Calculates the Vorticity with constant dx, dy, periodic BCs, and Ghost Cells
    
    Params:
    -------
    x, y    2D arrays of float, x and y direction points
    u, v    2D arrays of float, u and v velocity vectors
    
    Returns:
    -------
    Vorticity  2D array of float, curl of the velocity field
    """
    import numpy as np
    dV = np.empty_like(v, dtype=float)
    dU = np.empty_like(u, dtype=float)
    
    dV = np.gradient(v)
    dU = np.gradient(u)
    
    V_gradient = dV[1]/dx
    U_gradient = dU[0]/dy
    
    Vorticity = V_gradient - U_gradient
    Vorticity[-1,:] = Vorticity[0,:]
    Vorticity[:,-1] = Vorticity[:,0]
    Vorticity[-2,:] = Vorticity[-1,:]
    Vorticity[:,-2] = Vorticity[:,-1]
    
    return Vorticity