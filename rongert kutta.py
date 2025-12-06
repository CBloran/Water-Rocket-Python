stepsNbr = 300 # number of steps in the simulation

# Initial conditions
t0 = 0.00001 # Initial time
h0 = 0.00001 # Initial height
v0 = 0 # Initial velocity
Vw0 = 0.0005

def Rungekutta(t0, h0, v0, Vw0, stepNbr):
    
    def systeme_complet(t, y):
        """
        y = [h, v, Vw]
        Retourne [dh_dt, dv_dt, dVw_dt]
        """
        h, v, Vw = y
        
        # Équation de position
        dh_dt = v
        
        # Équation de volume d'eau
        if Vw > 0:  # Éviter les valeurs négatives
            dVw_dt = equation_water_volume(v, Vw)
            dv_dt = equation_vel
        else:
            dVw_dt = 0
            Vw = 0
        
        # Équation de vitesse
        dv_dt = equation_vel(v, Vw)

        
        
        return [dh_dt, dv_dt, dVw_dt]
    
    # Initialisation
    ts = [t0]
    ys = [[h0, v0, Vw0]]  # Stocker toutes les variables dans une liste
    
    deltaT = 0.01  # Pas de temps constant
    
    for i in range(stepNbr):
        t_current = ts[-1]
        y_current = ys[-1]
        print(y_current[1])
        
        # RK4 standard
        k1 = systeme_complet(t_current, y_current)
        
        k2 = systeme_complet(t_current + deltaT/2, 
                           [y_current[j] + deltaT/2 * k1[j] for j in range(3)])
        
        k3 = systeme_complet(t_current + deltaT/2,
                           [y_current[j] + deltaT/2 * k2[j] for j in range(3)])
        
        k4 = systeme_complet(t_current + deltaT,
                           [y_current[j] + deltaT * k3[j] for j in range(3)])
        
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) 
                for j in range(3)]
        
        # Conditions d'arrêt
        if y_new[0] <= 0:  # Au sol
            y_new[0] = 0
            y_new[1] = 0
        
            
        ts.append(t_current + deltaT)
        ys.append(y_new)
    
    # Séparation des résultats
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    Vws = [y[2] for y in ys]
    
    return ts, hs, vs, Vws