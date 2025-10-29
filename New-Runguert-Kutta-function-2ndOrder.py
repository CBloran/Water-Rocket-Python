def RungeKuttaSecondOrderGeneral(f, _Tfinal, _stepsNbr, _t0, _y0, _dy0=None):
    
    #======PARAMETERS=======#
    # f differential equation fonction
    # _Tfinal final time
    # _stepsNbr number of steps  
    #_t0 Initial time
    #_y0 Initial pos
    #_dy0 Initial speed leave blank if first order equation

    ts = [_t0]    
    ys = [_y0]    
    dys = [_dy0]  
    
    deltaT = _Tfinal / _stepsNbr
    
    for i in range(_stepsNbr):
        
        if _dy0 != None:
            m1_y = dys[-1]
            m1_dy = f(ts[-1], ys[-1], dys[-1])
        
            m2_y = dys[-1] + m1_dy * deltaT/2
            m2_dy = f(ts[-1] + deltaT/2, ys[-1] + m1_y * deltaT/2, dys[-1] + m1_dy * deltaT/2)
        
            m3_y = dys[-1] + m2_dy * deltaT/2
            m3_dy = f(ts[-1] + deltaT/2, ys[-1] + m2_y * deltaT/2, dys[-1] + m2_dy * deltaT/2)
        
            m4_y = dys[-1] + m3_dy * deltaT
            m4_dy = f(ts[-1] + deltaT, ys[-1] + m3_y * deltaT, dys[-1] + m3_dy * deltaT)

            m_y = (m1_y + 2*m2_y + 2*m3_y + m4_y) / 6
            m_dy = (m1_dy + 2*m2_dy + 2*m3_dy + m4_dy) / 6
        else:
            m1_y = f(ts[-1], ys[-1])
            m2_y = f(ts[-1] + deltaT / 2, ys[-1] + m1_y * deltaT / 2)
            m3_y = f(ts[-1] + deltaT / 2, ys[-1] + m2_y * deltaT / 2)
            m4_y = f(ts[-1] + deltaT, ys[-1] + m3_y * deltaT)
            m_y = (m1_y + 2 * m2_y + 2 * m3_y + m4_y) / 6

        
        
        
        next_y = ys[-1] + m_y * deltaT
        if _dy0 != None:
            next_dy = dys[-1] + m_dy * deltaT
            dys.append(next_dy)
        
        ts.append(ts[-1] + deltaT)
        ys.append(next_y)
        
    if _dy0 != None:
        return ts, ys, dys
    else:
        return ts, ys
    



def f(t, y, dy_dt=None):
    """ d2y/d2t + a*dy/dt + b*y = 0 """

    if dy_dt == None:
        dydt = -0.5*y
        return dydt
    else:
        dy2dt2 = -0.1 * dy_dt - 2**2 * y
        return dy2dt2
def g(t, y):

    dydt = -9.81*t + 20
    return dydt

# Utilisation
#ts, ys, dys = RungeKuttaSecondOrderGeneral(f, 2, 20, 0, 1, 0)
ts, ys = RungeKuttaSecondOrderGeneral(g, 4, 20, 0, 0)
print("Temps : ",ts)
print("Solutions y(t) : ",ys)
#print("Dérivées dy/dt : ",dys)