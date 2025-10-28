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
De = float(0.02)            # Diamètre de la buse (m)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
V = float(0.002)          # Volume total du rocket (m³)
p_ino = float(500000)           # Pression initiale à l'intérieur (Pa)
mb = float(1.7)             # Masse structurelle (kg)
Vwo = float(1.5)              # Volume d'eau initial dans la fusée
mw = rho_w * Vwo
k = Vwo /V

#retranscription des équations du document que j'ai (Alexandre) présenté, les arguments des fonctions
#permettent de connaître les paramètres / variables à mesurer/calculer

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
def FD_W(v):
    """
    v : la vitesse de la fusée
    Équations (4) :
    F_D = 1/2 * rho_atm * v² * Cd * A   -> force de traînée
    W   = (mb + mw) * g                 -> poids total du rocket
    """
    FD = 0.5 * rho_atm * v**2 * Cd * A
    W = (mb + mw) * g
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




#=======================================================
########           Graphique txt             ##########
#=======================================================

# fonction d'affichage des de liste de données suivant l'axe vertical parce que c est plus facil et que ca marche tout aussi bien


def txtGraph(xs: list, ys: list):
    source_file = open("Output.txt", "w")
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

Tfinal = 4 # lenght of the simulation
stepsNbr = 20 # number of steps in the simulation

# Initial conditions
y0 = 0 # Initial pos
t0 = 0 # Initial time

# differential equation dy/dt = f(t, y)
def f(t, y): 
    dydt = -9.81*t + 20

    return dydt

##### FUNCTION #####
def RungeKutta(_y0, _t0, _stepsNbr, _Tfinal):
    ts = [_t0] # Time list
    ys = [_y0] # Solution list

    deltaT = _Tfinal / _stepsNbr

    for i in range(_stepsNbr):

        m1 = f(ts[-1], ys[-1])
        m2 = f(ts[-1] + deltaT / 2, ys[-1] + m1 * deltaT / 2)
        m3 = f(ts[-1] + deltaT / 2, ys[-1] + m2 * deltaT / 2)
        m4 = f(ts[-1] + deltaT, ys[-1] + m3 * deltaT)
        m = (m1 + 2 * m2 + 2 * m3 + m4) / 6
        
        next_y = ys[-1] + m * deltaT

        ts.append(ts[-1] + deltaT)
        ys.append(next_y)


        if abs(ys[-2] - ys[-1]) < 0.001 or ys[-1] <= 0 :  # stop when the value change is negligible
            ys.pop(-1)
            ts.pop(-1)
            break

    print("ronguekutta", ts)
    print("ronguekutta", ys)
    x_simple = ts ##[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y_simple = ys ##[2, 4, 6, 8, 10, 8, 6, 4, 2, 0]
    
    # Appel avec couleur bleue
    ##graphique_temporel(x_simple, y_simple, 40, 40, 'bleu')
    ##graphique_temporel(x_simple, y_simple, 20, 20, 'vert')
    txtGraph(ts, ys)
    return (ts, ys)

#txtGraph()
tsys = RungeKutta(y0, t0, stepsNbr, Tfinal)
