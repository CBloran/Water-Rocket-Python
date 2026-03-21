
import math

mb = 0.308
rho_eau = 1000
P_atm = 101325
Cd = 0.71
rho_air = 1.225
D_tuyere = 0.008
D_fusee = 0.09
A_tuyere = math.pi * (D_tuyere**2) / 4
A_frontale =  math.pi * (D_fusee**2) / 4
V = 0.0015
g = 9.81


p0 = 3e5
mw0 = 0.0007



def vitesseEjection(P):
 
    delta_P = max(P - P_atm, 0.0)
 
    return math.sqrt(2.0 * delta_P / (rho_eau*(1-(A_tuyere/A_frontale)**2)))
 

def calculateF_grav(mw):

    return (mb + mw) * g

def calculateDrag(v):

    return 0.5 * rho_air * Cd * A_frontale * v **2

def Pression(mw):
    
    P = (p0*((V - mw0/rho_eau)**1.4))/((V - mw/rho_eau)**1.4)
    return P




def Euler(mw0, stepNbr=50000):


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

        F_poussee = 2*A_tuyere*(P - P_atm)

        dv_dt = (F_poussee - F_grav - F_trainee)/(mb + mw) 

        dmw_dt = -rho_eau * A_tuyere * v_e

        return [dh_dt, dv_dt, dmw_dt], (F_poussee - F_grav)


    while ys[-1][2] > 0:
        
        t_current = ts[-1]
        y_current = ys[-1]
        print(y_current)
        CalculatedVal = []

        if y_current[2] > 0.00000001:
                
            k, F = systemeCompletEau(y_current[:3])
            CalculatedVal.append(F)
            y_new = [y_current[i] + deltaT * k[i] for i in range(3)] + CalculatedVal

            ys.append(y_new)
            ts.append(t_current + deltaT)
        
        elif y_current[0] > 0.00000001:
            k, F = systemeCompletEau(y_current[:3])
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
        