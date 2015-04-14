def F1V(u_P, u_E, u_W, u_N, u_S, dx, dy, nu):
    """Computes the Viscous terms
    
    Params:
    -------
    u_P, u_E, u_W, u_N, u_S   float, u velocity points
    dx, dy                    float, x and y step size
    nu                        float, fluid viscosity
    
    Returns:
    -------
    F1_V                      float, u direction viscous terms
    """
    
    F1_V = ((nu/dx)*(((u_E - u_P)/dx) - ((u_P - u_W)/dx))) + ((nu/dy)*(((u_N -\
        u_P)/dy) - ((u_P - u_S)/dy)))
    
    return F1_V

def F2V(v_P, v_E, v_W, v_N, v_S, dx, dy, nu):
    """Computes the Viscous terms
    
    Params:
    -------
    v_P, v_E, v_W, v_N, v_S   float, v velocity points
    dx, dy                    float, x and y step size
    nu                        float, fluid viscosity
    
    Returns:
    -------
    F2_V                      float, v direction viscous terms
    """
    
    F2_V = ((nu/dx)*(((v_E - v_P)/dx) - ((v_P - v_W)/dx))) + ((nu/dy)*(((v_N -\
        v_P)/dy) - ((v_P - v_S)/dy)))
    
    return F2_V