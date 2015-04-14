def FP1(P_e, P_w, dx):
    """Computes the Pressure terms
    
    Params:
    -------
    P_e, P_w                float, Pressure Points
    dx, dy                  float, x and y step size
    
    Returns:
    --------
    F1_P                    float, u direction Pressure terms
    """
    
    F1_P = (P_e + P_w)/dx
    
    return F1_P

def FP2(P_n, P_s, dy):
    """Computes the Pressure terms
    
    Params:
    -------
    P_n, P_s                float, Pressure Points
    dx, dy                  float, x and y step size
    
    Returns:
    --------
    F2_P                    float, v direction Pressure terms
    """
    
    F2_P = (P_n + P_s)/dy
    
    return F2_P