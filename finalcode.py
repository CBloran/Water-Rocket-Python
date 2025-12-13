
import math

# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.81               # Accélération de la gravité (m/s²)
rho_w = 1.0e3         # Densité de l'eau (kg/m³)
rho_atm = 1.23        # Densité de l'air (kg/m³)
Cd = 0.35             # Coefficient de traînée aérodynamique
gamma = 1.4           # Coefficient adiabatique de l'air
patm = 101325         # Pression atmosphérique (Pa)

# Paramètres géométriques et initiaux du rocket
#les noms avec o sont des paramètres à t = 0s (on peut pas mettre 0 dans la variable donc o = initial, tandis que "in" comme dans p_ino est pour "intérieur" ou "interne")
# les paramètres qu'on ne connait pas sont notés float pour l'instant
D = float(0.1)               # Diamètre de la fusée (m)
De = float(0.008)            # Diamètre de la buse (m)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
V = float(0.0015)          # Volume total du rocket (m³)
p_ino = float(300000)           # Pression initiale à l'intérieur (Pa)
#p_ino += 101325             # Convertir la pression relative en pression absolue
mb = float(0.1)             # Masse structurelle (kg)
Vwo = float(0.0005)              # Volume d'eau initial dans la fusée
mw = rho_w * Vwo
k = Vwo /V
k0 = Vwo/V


# ============================================================
#  FORCES EN PRESENCE
# ============================================================
def Weight(mw):
    """
    mw : water mass
    Calculate the weight of the rocket
    """
    
    W = (mb + mw) * g
    return W
def Drag(v):
    """
    v : Rocket speed
    Calcculate the drag force
    """
    FD = 0.5 * rho_atm * v**2 * Cd * A
    return FD
def Thrust(v_e):
    """
    v_e : exit velocity of water
    calculate the thrust force based on the exit velocity of water
    """
    F_thrust = rho_w * Ae * v_e**2
    return F_thrust
# ============================================================
def internal_pressure(Vw):
    """
    Vw : volume d'eau restant
    Calcul de la pression interne en fonction du volume d'eau restant.
    """
    Va = V - Vw  # remaining air volume
    Vao = V - Vwo  # initial air volume

    if Va > 0:
        p_in = p_ino * (Vao / Va)**gamma # calculate the remaining internal pressure based on adiabatic law
        return p_in
    else:
        return patm 

def water_exit_velocity(k, p_in):
    """
    k : ratio of remaining water volume to total volume
    p_in : internal pressure
    """
    try:
        v_e = ((2*(p_ino*((1-k0)/(1-k))**(gamma)-patm))/((rho_w)*(1-((Ae)/(A))**2)))**(1/2) #calculate the exit velocity of water based on bernouilli's equation
        return v_e
    except:
        return 0

def equation_vel(v, Vw):
            
            mw_current = rho_w * Vw  # Masse d'eau actuelle
            F_drag = Drag(v)
            F_weight = Weight(mw_current)  
            
            if Vw > 0:
                p_in = internal_pressure(Vw)
                v_e = water_exit_velocity(Vw/V, p_in)
                F_thrust = Thrust(v_e)
                #print(F_thrust)
            else:
                F_thrust = 0
            
            if Vw > 0:
                dv_dt = (F_thrust - F_drag - F_weight) / (mb + mw_current)
            else:
                dv_dt = (-F_drag - F_weight) / mb
            return dv_dt

def equation_water_volume(v, Vw):

    k = Vw / V                            # calculate the ratio of remaining water 
    p_in = internal_pressure(Vw)
    v_e = water_exit_velocity(k, p_in)
    dVw_dt = -Ae * v_e        # rate of change of water volume
    return dVw_dt

#=======================================================
########           Graphique txt             ##########
#=======================================================

# fonction d'affichage des de liste de données suivant l'axe vertical parce que c est plus facil et que ca marche tout aussi bien


def txtGraph(xs: list, ys: list, FileName = "Output.txt"):
    source_file = open(FileName, "w")
    x_simple = xs[:] ##[1, 2.1, 3.6, 4.4, 5.5, 6, 7, 8, 9, 10]  test data
    y_simple = ys[:] ##[2.1, 4, 6, 8, 10, 8.5, 6.6, 4.8, 2.8, 0] test data

    for i in range(len(x_simple)):       # only keep the rounded value of our data because we cant use decimal number to set the position of the text
        x_simple[i] = round(x_simple[i])
        y_simple[i] = round(y_simple[i])

    max_x = max(x_simple)                
    max_y = max(y_simple)                # find the max value of the list to adjust the pos of the text    
    

    for i in range(len(y_simple)):
        line = (y_simple[i] - 1) * " " + str(round(ys[i], 1)) + (max_y - y_simple[i]) * " "
        source_file.write(line + "\n")



#=======================================================
########           Runguert-Kutta             ##########
#=======================================================
##### Simulation parameters #####




def Rungekutta(stepNbr=3000):
    """
    stepNbr : number of points in the simulation
    """

    # Initial conditions
    t0 = 0.00001 # Initial time
    h0 = 0.00001 # Initial height
    v0 = 0 # Initial velocity
    Vw0 = 0.0005
    
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
    i = 0
    
    while ys[-1][0] > 0 and i < stepNbr:
        i += 1
        t_current = ts[-1]
        y_current = ys[-1]
        print(y_current[1])
        
        # RK4 standard
        k1 = systeme_complet(t_current, y_current)
        
        k2 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k1[j] for j in range(3)])
        
        k3 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k2[j] for j in range(3)])
        
        k4 = systeme_complet(t_current + deltaT, [y_current[j] + deltaT * k3[j] for j in range(3)])
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(3)]
            
        ts.append(t_current + deltaT)
        ys.append(y_new)
    
    ys.append([0, 0, 0])
    
    # Séparation des résultats
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    Vws = [y[2] for y in ys]
    
    return ts, hs, vs, Vws





tsys = Rungekutta()
tsys = Rungekutta()
#print(tsys[2])
txtGraph(tsys[0], tsys[1], "Height.txt")