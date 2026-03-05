
import math
import turtle
import matplotlib.pyplot as plt

# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.81               # Accélération de la gravité (m/s²)
rho_w = 1.0e3         # Densité de l'eau (kg/m³)
rho_atm = 1.23        # Densité de l'air (kg/m³)
Cd = 0.73             # Coefficient de traînée aérodynamique
gamma = 1.4           # Coefficient adiabatique de l'air
patm = 101325         # Pression atmosphérique (Pa)
R_air = 287.05

# Paramètres géométriques et initiaux du rocket
#les noms avec o sont des paramètres à t = 0s (on peut pas mettre 0 dans la variable donc o = initial, tandis que "in" comme dans p_ino est pour "intérieur" ou "interne")
# les paramètres qu'on ne connait pas sont notés float pour l'instant
D = float(0.1)               # Diamètre de la fusée (m)
De = float(0.008)            # Diamètre de la buse (m)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
V = float(0.0015)          # Volume total du rocket (m³)
p_ino = float(360000)           # Pression initiale à l'intérieur (Pa)
T_initial = 293.15  # 20°C
#p_ino += 100000             # Convertir la pression relative en pression absolue
mb = float(0.3)             # Masse structurelle (kg)
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
        return 0 

def water_exit_velocity(k, p_in):
    """
    k : ratio of remaining water volume to total volume
    p_in : internal pressure
    """
    try:
        v_e = ((2*(p_ino*((1-k0)/(1-k))**(gamma))-patm)/((rho_w)*(1-((Ae)/(A))**2)))**(1/2) #calculate the exit velocity of water based on bernouilli's equation
        #  OTHER BERNOUILLI EQUATION  #
        #delta_P = p_in - patm
        #v_e = (2.0 * delta_P / (rho_w*(1-(Ae/A)**2)))**(1/2) #calculate the exit velocity of water based on bernouilli's equation
        
        return v_e
    except:
        return 0

def equation_vel(v, Vw, p_in):
            
            mw_current = rho_w * Vw  # Masse d'eau actuelle
            F_drag = Drag(v)
            F_weight = Weight(mw_current)  
            
            if Vw > 0:
                
                v_e = water_exit_velocity(Vw/V, p_in)
                F_thrust = Thrust(v_e)
                #print(F_thrust)
            else:
                F_thrust = 0
            
            if Vw > 0:
                dv_dt = (F_thrust - F_drag - F_weight) / (mb + mw_current)
            else:
                dv_dt = (-F_drag - F_weight) / mb
            return dv_dt, (F_thrust)
            
def equation_vel_air(v, Vw, p_in):
            
            F_drag = Drag(v)
            F_weight = Weight(0)
            T_air = T_initial * (p_in / p_ino) ** ((gamma - 1) / gamma)
            
            if p_in > patm:
                pressure_ratio = p_in / patm
                
                # Vitesse d'éjection
                v_e = math.sqrt(2 * gamma / (gamma - 1) * R_air * T_air * (1 - (1 / pressure_ratio) ** ((gamma - 1) / gamma)))
                F_thrust = Ae * v_e**2 
                
            else:
                F_thrust = 0
            
            dv_dt = (0 - F_drag - F_weight) / mb
            
            return dv_dt, F_thrust

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
########           Graph V2 MathPlot             #########
#=======================================================
def graphMathPlot(Times, dataIn, additionalValues=None):
    

    data = dataIn

    x = Times 

    plt.figure(figsize=(10, 5))

    # Séparation des composantes
    h  = [row[0] for row in data]
    plt.plot(x, h,  label="h",  color="red",   marker="o", markersize = 2)

    if len(data[0]) > 1:
        v  = [row[1] for row in data]
        plt.plot(x, v,  label="v",  color="blue",  marker="o", markersize = 2)

    if len(data[0]) > 2:
        Vw = [row[2]*1000 for row in data]
        plt.plot(x, Vw, label="Vw", color="green", marker="o", markersize = 2)

    if len(data[0]) > 3:
        p_in = [row[3]/1000 for row in data]
        plt.plot(x, p_in, label="p_in", color="purple", marker="o", markersize = 2)
    
    if len(data[0]) > 4:
        F = [row[4] for row in data]
        plt.plot(x, F, label="F", color="orange", marker="o", markersize = 2)
    

    
    
    
    
    if additionalValues:
        plt.plot(x, additionalValues, label="F_thrust", color="orange", marker="o", markersize = 2)

    plt.xlabel("Index")
    plt.ylabel("Valeur")
    plt.title("Évolution de h, v et Vw")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()



#=======================================================
########           Runguert-Kutta             ##########
#=======================================================
##### Simulation parameters #####




def Rungekutta(stepNbr=5000):
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
        Params:
        t = time
        y = [h, v, Vw, p_in, F]
        return [dh_dt, dv_dt, dVw_dt, dp_in_dt, F]
        """
        h, v, Vw, p_in, F = y
        
        # Équation de position
        dh_dt = v
        
        # Équation de volume d'eau
        if Vw > 0.0000000001:  # Éviter les valeurs négatives
            dVw_dt = equation_water_volume(v, Vw)
            
        else:
            dVw_dt = 0
            Vw = 0
        
        # Équation de vitesse
        if Vw > 0.0000000001:
            p_in_next = internal_pressure(Vw)
            dp_in_dt = (p_in_next - p_in) / deltaT
            dv_dt, F = equation_vel(v, Vw, p_in)
            
        else:
            dh_dt = v
            if p_in > patm: 
                dp_in_dt = -70000*(1/t)
            else:
                 dp_in_dt = 0
                 p_in = patm
            
            dv_dt, F = equation_vel_air(v, Vw, p_in)

        return [dh_dt, dv_dt, dVw_dt, dp_in_dt, F]
    
    
    # Initialisation
    ts = [t0]
    Fs = [0]
    ys = [[h0, v0, Vw0, p_ino, 0]]  # Stocker toutes les variables dans une liste
    
    deltaT = 0.01  # Pas de temps constant
    i = 0
    
    while ys[-1][0] > 0 and i < stepNbr :
        i += 1
        t_current = ts[-1]
        y_current = ys[-1]
        #print(y_current[0])
        
        p_in = internal_pressure(y_current[2])
        v_e = water_exit_velocity(y_current[2]/V, p_in)
        F_thrust = Thrust(v_e)
        Fs.append(p_in)
        #print(p_in)
        # RK4 standard
        k1 = systeme_complet(t_current, y_current)
        
        k2 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k1[j] for j in range(4)]+[y_current[4]])
        
        k3 = systeme_complet(t_current + deltaT/2, [y_current[j] + deltaT/2 * k2[j] for j in range(4)]+[y_current[4]])
        
        k4 = systeme_complet(t_current + deltaT, [y_current[j] + deltaT * k3[j] for j in range(4)]+[y_current[4]])
        
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(4)]+[k1[4]]
            
        ts.append(t_current + deltaT)
        ys.append(y_new)


    
    
    graphMathPlot(ts, ys)
    # Séparation des résultats
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    Vws = [y[2] for y in ys]
    Fs = [y[4] for y in ys]
    print(Fs)
    return ts, hs, vs, Vws





tsys = Rungekutta()

#print(tsys[2])
txtGraph(tsys[0], tsys[1], "Height.txt")