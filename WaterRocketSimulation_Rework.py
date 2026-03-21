
import math

mb = 0.308
rho_eau = 1000
P_atm = 101325
Cd = 0.71
rho_air = 1.225
D_tuyere = 0.009
D_fusee = 0.09
A_tuyere = math.pi * (D_tuyere**2) / 4
A_frontale =  math.pi * (D_fusee**2) / 4
V = 0.0015
g = 9.81


p0 = 3e5 + P_atm
mw0 = 0.7



def vitesseEjection(P):
 
    delta_P = max(P - P_atm, 0.0)
 
    return math.sqrt(2.0 * delta_P / (rho_eau*(1-(A_tuyere/A_frontale)**2)))
 

def calculateF_grav(mw):
    F_grav = (mb + mw) * g
    if mw > 0.001:
        print(F_grav)
    return F_grav

def calculateDrag(v):

    return 0.5 * rho_air * Cd * A_frontale * v **2

def Pression(mw):
    
        #P = (p0*((V - mw0/rho_eau)**1.4))/((V - mw/rho_eau)**1.4)
        P = p0 * ((V - mw0/rho_eau)**1.4)/((V - mw/rho_eau)**1.4)
        return P




def Euler(mw0, stepNbr=500000):


    h0 = 0
    v0 = 0
    ts = [0]

    ys = [[h0, v0, mw0, 0]]
    ts = [0]
    deltaT = 0.001

    
    def systemeCompletEau(y):
        
        h, v, mw, = y

        dh_dt = v
        
        P = Pression(mw)
        F_grav = calculateF_grav(mw)
        F_trainee = calculateDrag(v)
        v_e = vitesseEjection(P)

        dmw_dt = -rho_eau * A_tuyere * v_e

        F_poussee = -dmw_dt*v_e

        dv_dt = (F_poussee - F_grav - F_trainee)/(mb + mw) 

        

        return [dh_dt, dv_dt, dmw_dt], (F_poussee - F_grav)
    
    def systemeCompletAutre(y):
        
        h, v, mw, = y

        dh_dt = v
        
        P = P_atm
        F_grav = calculateF_grav(mw)
        F_trainee = calculateDrag(v)
        v_e = 0

        F_poussee = 0

        dv_dt = (F_poussee - F_grav - F_trainee)/(mb) 

        dmw_dt = 0

        return [dh_dt, dv_dt, dmw_dt], (-F_grav)


    while ys[-1][0] >= 0:
        
        t_current = ts[-1]
        y_current = ys[-1]
        
        CalculatedVal = []

        if y_current[2] > 0.00000001:
                
            k, F = systemeCompletEau(y_current[:3])
            CalculatedVal.append(F)
            y_new = [y_current[i] + deltaT * k[i] for i in range(3)] + CalculatedVal
            #print(y_current)
            ys.append(y_new)
            ts.append(t_current + deltaT)
        
        else:
            k, F = systemeCompletAutre(y_current[:3])
            CalculatedVal.append(F)
            y_new = [y_current[i] + deltaT * k[i] for i in range(3)] + CalculatedVal
            
            ys.append(y_new)
            ts.append(t_current + deltaT)
    
    with open("resultats.txt", "w", encoding="utf-8")as w:
        for i in range(len(ys)):
            line =  str(ts[i])
            for val in ys[i]:
                line += "," + str(val)
            line += "\n"
            w.write(line)
    

Euler(0.7)
        