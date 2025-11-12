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
p_ino = float(500000)           # Pression initiale à l'intérieur (Pa)
mb = float(0.1)             # Masse structurelle (kg)
Vwo = float(0.0005)              # Volume d'eau initial dans la fusée
mw = rho_w * Vwo
k = Vwo /V
k0 = Vwo/V

#retranscription des équations du document que j'ai (Alexandre) présenté, les arguments des fonctions
#permettent de connaître les paramètres / variables à mesurer/calculer
v_e=((2*(p_ino*((1-k0)/(1-k))**(gamma)-patm))/((rho_w)*(1-((Ae)/(A))**2)))**(1/2)
# ============================================================
# (1) ÉQUATION DE VARIATION DE LA MASSE D'EAU
# ============================================================
def dmw_dt(v_e):
    """
    v_e : vitesse d'éjecion de l'eau dans la buse 
    Équation (1) du papier :
    dm_w/dt = -rho_w * Ae * v_e
    -> décrit la variation de la masse d’eau à mesure qu’elle est expulsée.
    """
    return -rho_w * Ae * v_e

# ============================================================
# (3) MASSE D’EAU RESTANTE DANS LE ROCKET
# ============================================================
def mw_remaining(delta_mw):
    """
    delta_mw : le delta de masse d'eau
    Équation (3) :
    m_w = rho_w * (k * V0 - ΔV)
    -> masse d’eau restante à l’instant t.
    """
    return rho_w * k * V - delta_mw


# ============================================================
# (4) FORCES DE TRAÎNÉE ET DE GRAVITÉ
# ============================================================
def FD_W(v, mw1 = mw):
    """
    v : la vitesse de la fusée
    Équations (4) :
    F_D = 1/2 * rho_atm * v² * Cd * A   -> force de traînée
    W   = (mb + mw) * g                 -> poids total du rocket
    """
    FD = 0.5 * rho_atm * v**2 * Cd * A
    W = (mb + mw1) * g
    return FD, W


# ============================================================
# (5) MASSE ET PRESSION DU GAZ INTERNE
# ============================================================
def gazi(Va):
    """
    Va : le volume d'air à un instant t
    Équation (5) :
    p_in = p_ino * (Vao / Va)^gamma
    -> décrit la diminution de la pression interne selon la loi adiabatique.
    Vao = V - Vw
    -> le volume d'air initial vaut le total  moins le volume d'eau
    rho_ino = (rho_w*Vw+rho_atm*Va)/(Vw+Va)
    -> calcul du rho initial du mélange eau-air dans le réservoir
    ma = rho_ino * (1-k)**V
    -> calcul de la masse d'air
    """
    Vao = V - Vwo # Volume d'air restant initial

    rho_ino = (rho_w * Vwo + rho_atm * Vao) / (Vwo + Vao)
    return p_ino * (Vao / Va)**gamma, rho_ino * (1-k)**Vwo

# ============================================================
# (6) SYSTÈME D'ÉQUATIONS DIFFÉRENTIELLES DE LA PHASE PROPULSIVE
# ============================================================
def rocket_dynamics(v, v_e):
    """
    mêmes définitions de v et v_e que précédemment
    Équations (6) :
        dm_w/dt   = -rho_w * Ae * v_e
        dv/dt     = (F - F_D - W) / (mb + mw)
        dp_in/dt  = -γ * p_in * v_e * Ae / (k * V0)
    F = -dmdt*v_e
    """

    # --- Dérivée de la masse d’eau
    dmdt = -rho_w * Ae * v_e

    # --- Poussée instantanée (force exercée par l’eau sur le rocket)
    #F = -dmdt*v_e //old
    F = rho_w * Ae * v_e**2 # force de poussée en fonction de v_e plus interescent pour la resolution de l'equation differentielle

    # --- Forces aérodynamiques
    FD, W = FD_W(v)

    # --- Accélération du rocket
    dvdt = (F - FD - W) / (mb + mw)

    # --- Variation de la pression interne
    dp_indt = -gamma * p_ino**((1+gamma)/gamma) * v_e * Ae / ((1-k) * V * p_ino**((1+gamma)/gamma))

    return F, dvdt, dp_indt

########################################################################################################
######################    WORKING ON A COMPLETE REWORK OF THE DIFFERENTIAL EQUATIONS SYSTEM   ##########
########################################################################################################
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
        #v_e = math.sqrt((2 * (p_in - patm)) / (rho_w * (1 - (Ae / A)**4))) #calculate the exit velocity of water based on bernouilli's equation
        #v_e = math.sqrt((2 * (patm - p_in)) / (rho_w * ((Ae/A)**4 - 1))) #calculate the exit velocity of water based on bernouilli's equation
        v_e=((2*(p_ino*((1-k0)/(1-k))**(gamma)-patm))/((rho_w)*(1-((Ae)/(A))**2)))**(1/2)
        return v_e
    except:
        return 0
    


def thrust_phase_system(t, v, Vw):
    """
    t : time
    v : speed
    Vw : remaining water volume
    differential equations system during the thrust phase
    """

    if Vw <= 0 and v <= 0:
        return 0, 0, 0, 0  # Fin de la phase propulsive lorsque tout l'eau est expulsée
    
    p_in = internal_pressure(Vw)          # calculate the internal pressure
    k = Vw / V                            # calculate the ratio of remaining water volume to total volume
    v_e = water_exit_velocity(k, p_in)    # calculate the exit velocity of water

    mw = rho_w * Vw                       # calculate the remaining water mass
    if Vw > 0:
        F_thrust = rho_w * Ae * v_e**2        # calculate the thrust force
    else:
        F_thrust = 0

    F_drag, F_weight = FD_W(v)            # calculate drag force and weight

    if Vw > 0:
        dh_dt = v                 # rate of change of height
        dv_dt = (F_thrust - F_drag - F_weight) / (mb + mw)  # rate of change of velocity
        dVw_dt = -Ae * v_e        # rate of change of water volume
    else:
        dh_dt = v
        dv_dt = ( - F_drag - F_weight) / mb  # rate of change of velocity without thrust
        dVw_dt = 0                 # no more water to expel
    return dh_dt, dv_dt, dVw_dt

#=========================================================================
########           Équations différentielles separer            ##########
#=========================================================================
def equation_pos(V, Vw):

    dh_dt = V
    return dh_dt

def equation_vel(v, Vw):
            
            mw_current = rho_w * Vw  # Masse d'eau actuelle
            F_drag, F_weight = FD_W(v, mw_current)  # Utiliser masse actuelle
        
            if Vw > 0:
                p_in = internal_pressure(Vw)
                v_e = water_exit_velocity(Vw/V, p_in)
                F_thrust = rho_w * Ae * v_e**2
                print(F_thrust)
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


stepsNbr = 300 # number of steps in the simulation

# Initial conditions
t0 = 0.00001 # Initial time
h0 = 0.00001 # Initial height
v0 = 0 # Initial velocity
Vw0 = 0.0005

# differential equation dy/dt = f(t, y)
def f(t, y): 
    dydt = -9.81*t + 20

    return dydt


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
            
        if y_new[2] <= 0:  # Plus d'eau
            y_new[2] = 0
            
        ts.append(t_current + deltaT)
        ys.append(y_new)
    
    # Séparation des résultats
    hs = [y[0] for y in ys]
    vs = [y[1] for y in ys]
    Vws = [y[2] for y in ys]
    
    return ts, hs, vs, Vws





tsys = Rungekutta(t0, h0, v0, Vw0, stepsNbr)
#print(tsys[2])
txtGraph(tsys[0], tsys[1], "Height.txt")