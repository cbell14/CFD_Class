def Pressure_Terms(P, dx, dy):
    """Computes the Pressure terms
    
    Params:
    -------
    P                       float, Pressure Points
    dx, dy                  float, x and y step size
    
    Returns:
    --------
    F1_P                    float, u direction Pressure terms
    F2_P                    float, v direction Pressure terms
    """
    
    F1_P = (P[:,1:] - P[:,:-1])/dx
    F2_P = (P[1:,:] - P[:-1,:])/dy
    
    return F1_P, F2_P
