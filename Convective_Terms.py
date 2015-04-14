def F1C(u_e, u_w, u_n, u_s, v_n, v_s, dx, dy):
    """Computes the Convec terms
    
    Params:
    -------
    u_e, u_w, u_n, u_s        float, u velocity points
    v_n, v_s                  float, v velocity points
    dx, dy                    float, x and y step size
    
    Returns:
    -------
    F1_C                      float, u direction convective terms
    """

    F1_C = -(u_e**2 - u_w**2)/dx - ((u_n*v_n) - (u_s*v_s))/dy
    
    return F1_C

def F2C(u_e, u_w, v_e, v_w, v_n, v_s, dx, dy):
    """Computes the Convec terms
    
    Params:
    -------
    u_e, u_w                  float, u velocity points
    v_e, v_w, v_n, v_s        float, v velocity points
    dx, dy                    float, x and y step size
    
    Returns:
    -------
    F2_C                      float, v direction convective terms
    """

    F2_C = -((u_e*v_e) - (u_w*v_w))/dx - (v_n**2 - v_s**2)/dy
    
    return F2_C