import math

# ============================================================
# CONSTANTES PHYSIQUES (À AJOUTER)
# ============================================================
g = 9.81
rho_w = 1000
rho_atm = 1.23
Cd = 0.35      # coefficient de trainer
gamma = 1.4
patm = 101325

# Paramètres fusée (À AJOUTER)
D = 0.2                     # Diamètre de la fusée [m]
De = 0.02                   # Diamètre de la buse [m]
A = math.pi * (D/2)**2      # Aire frontale [m²]
Ae = math.pi * (De/2)**2    # Aire de la buse [m²]
V = 0.001                  # Volume total du réservoir [m³] (1.5 litres)
p_ino = 500000              # Pression initiale [Pa] (~5 bars)
mb = 0.1                    # Masse à vide de la fusée [kg] - À AJUSTER
Vwo = 0.00035                 # Volume d'eau initial [m³] (0.5 litre)

def FD_W(v):
    """Calcule traînée et poids"""
    FD = 0.5 * rho_atm * v**2 * Cd * A
    W = mb * g
    return FD, W

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
        v_e = math.sqrt((2 * (p_in - patm)) / (rho_w * (1 - (Ae / A)**2))) #calculate the exit velocity of water based on bernouilli's equation
        return v_e
    except:
        return 0

# =========================================================================
# ÉQUATIONS DIFFÉRENTIELLES SÉPARÉES (CORRIGÉES)
# =========================================================================
def equation_pos(v, Vw):
    """dh/dt = v"""
    return v

def equation_vel(v, Vw):
    """dv/dt = f(v, Vw)"""
    if Vw <= 0 and v <= 0:
        return 0
    
    p_in = internal_pressure(Vw)
    k = Vw / V
    v_e = water_exit_velocity(k, p_in)
    mw = rho_w * Vw

    if Vw > 0:
        F_thrust = rho_w * Ae * v_e**2
    else:
        F_thrust = 0

    F_drag, F_weight_base = FD_W(v)
    F_weight = F_weight_base + mw * g  # Ajouter le poids de l'eau

    if Vw > 0:
        return (F_thrust - F_drag - F_weight) / (mb + mw)
    else:
        return (- F_drag - F_weight) / mb

def equation_water_volume(v, Vw):
    """dVw/dt = f(v, Vw)"""
    if Vw <= 0:
        return 0
    
    p_in = internal_pressure(Vw)
    k = Vw / V
    v_e = water_exit_velocity(k, p_in)

    if Vw > 0:
        return -Ae * v_e
    else:
        return 0
    
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

# =========================================================================
# RUNGE-KUTTA CORRIGÉ
# =========================================================================
# Version ultra-simplifiée si vous voulez garder votre structure actuelle
def Rungekutta_simple_until_landing(t0, h0, v0, Vw0):
    """
    Version très simple - simulation jusqu'au sol
    """
    def systeme_complet(h, v, Vw):
        return [
            equation_pos(v, Vw),
            equation_vel(v, Vw),  
            equation_water_volume(v, Vw)
        ]
    
    ts = [t0]
    hs = [h0]
    vs = [v0]
    Vws = [Vw0]
    
    deltaT = 0.01  # Pas de temps fixe
    max_steps = 100  # Sécurité

    for i in range(max_steps):
        t_current = ts[-1]
        h_current = hs[-1]
        v_current = vs[-1]
        Vw_current = Vws[-1]

        # Runge-Kutta
        k1 = systeme_complet(h_current, v_current, Vw_current)
        k2 = systeme_complet(
            h_current + deltaT/2 * k1[0], 
            v_current + deltaT/2 * k1[1], 
            Vw_current + deltaT/2 * k1[2]
        )
        k3 = systeme_complet(
            h_current + deltaT/2 * k2[0],
            v_current + deltaT/2 * k2[1], 
            Vw_current + deltaT/2 * k2[2]
        )
        k4 = systeme_complet(
            h_current + deltaT * k3[0],
            v_current + deltaT * k3[1],
            Vw_current + deltaT * k3[2]
        )

        m_h = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
        m_v = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
        m_Vw = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6

        new_h = h_current + deltaT * m_h
        new_v = v_current + deltaT * m_v
        new_Vw = Vw_current + deltaT * m_Vw

        # Arrêter si au sol
        if new_h <= 0:
            new_h = 0
            new_v = 0
            hs.append(new_h)
            vs.append(new_v)
            Vws.append(new_Vw)
            ts.append(t_current + deltaT)
            break

        hs.append(new_h)
        vs.append(new_v)
        Vws.append(new_Vw)
        ts.append(t_current + deltaT)

    return ts, hs, vs, Vws

# Utilisation :
ts, hs, vs, Vws = Rungekutta_simple_until_landing(0, 0.001, 0, Vwo)
txtGraph(ts, hs)