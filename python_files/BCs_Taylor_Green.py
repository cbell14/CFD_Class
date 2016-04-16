def BCs(u, v):
    """Applies Dirichlet & Neumann Boundary conditions on u and v from
    
    Params:
    ------
    u, v            2D arrays of float, calculated u and v
    
    Returns:
    -------
    u, v            2D array of float, u and v with BCs applied
    """
    u[:,0] = u[:,2]
    u[:,-1] = u[:,-3]
    u[0,:] = - u[1,:]
    u[-1,:] = - u[-2,:]

    v[:,0] = - v[:,1]
    v[:,-1] = - v[:,-2]
    v[0,:] = v[2,:]
    v[-1,:] = v[-3,:]

    return u, v