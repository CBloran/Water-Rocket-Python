import math

# ============================================================
# CONSTANTES PHYSIQUES (À AJOUTER)
# ============================================================
g = 9.81
rho_w = 1000
rho_atm = 1.23
Cd = 0.35
gamma = 1.4
patm = 101325

# Paramètres fusée (À AJOUTER)
D = 0.1
De = 0.02
A = math.pi * (D/2)**2
Ae = math.pi * (De/2)**2
V = 0.002
p_ino = 500000
mb = 1.7
Vwo = 1

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
def Rungekutta(t0, h0, v0, Vw0, stepNbr, Tfinal):
    
    def systeme_complet(h, v, Vw):
        """
        Retourne [dh_dt, dv_dt, dVw_dt]
        """
        return [
            equation_pos(v, Vw),           # dh/dt
            equation_vel(v, Vw),           # dv/dt  
            equation_water_volume(v, Vw)   # dVw/dt
        ]
    
    ts = [t0]
    hs = [h0]
    vs = [v0]
    Vws = [Vw0]

    deltaT = Tfinal / stepNbr
    
    print(f"DEBUG: Début simulation - deltaT={deltaT}")
    print(f"DEBUG: Conditions initiales - h0={h0}, v0={v0}, Vw0={Vw0}")

    for i in range(stepNbr):
        t_current = ts[-1]
        h_current = hs[-1]
        v_current = vs[-1]
        Vw_current = Vws[-1]
        
        print(f"DEBUG: Step {i} - t={t_current:.3f}, h={h_current:.3f}, v={v_current:.3f}, Vw={Vw_current:.6f}")

        # POINT 1
        k1 = systeme_complet(h_current, v_current, Vw_current)
        print(f"DEBUG: k1 = {k1}")
        
        # POINT 2
        k2 = systeme_complet(
            h_current + deltaT/2 * k1[0], 
            v_current + deltaT/2 * k1[1], 
            Vw_current + deltaT/2 * k1[2]
        )
        
        # POINT 3
        k3 = systeme_complet(
            h_current + deltaT/2 * k2[0],
            v_current + deltaT/2 * k2[1], 
            Vw_current + deltaT/2 * k2[2]
        )
        
        # POINT 4
        k4 = systeme_complet(
            h_current + deltaT * k3[0],
            v_current + deltaT * k3[1],
            Vw_current + deltaT * k3[2]
        )

        # MOYENNE PONDÉRÉE
        m_h = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
        m_v = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
        m_Vw = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6

        # MISE À JOUR
        new_h = h_current + deltaT * m_h
        new_v = v_current + deltaT * m_v
        new_Vw = Vw_current + deltaT * m_Vw
        
        print(f"DEBUG: Nouveaux - h={new_h:.3f}, v={new_v:.3f}, Vw={new_Vw:.6f}")

        hs.append(new_h)
        vs.append(new_v)
        Vws.append(new_Vw)
        ts.append(t_current + deltaT)

        # Conditions d'arrêt
        if new_h <= 0:  # Au sol
            print("DEBUG: Impact au sol - arrêt")
            hs[-1] = 0
            vs[-1] = 0
            break
            
        if new_Vw <= 0:  # Plus d'eau
            print("DEBUG: Plus d'eau - arrêt")
            Vws[-1] = 0
            break

    print(f"DEBUG: Simulation terminée - {len(ts)} points calculés")
    return ts, hs, vs, Vws

# Test avec des paramètres plus réalistes
Tfinal = 2  # Réduire le temps
stepsNbr = 50  # Réduire le nombre d'étapes pour voir ce qui se passe

# Conditions initiales
t0 = 0.0
h0 = 0.001    # Légèrement au dessus du sol
v0 = 0.0
Vw0 = 1   # Volume d'eau réaliste

print("=== DÉBUT SIMULATION ===")
ts, hs, vs, Vws = Rungekutta(t0, h0, v0, Vw0, stepsNbr, Tfinal)

if len(ts) > 1:
    print(f"\n=== RÉSULTATS ===")
    print(f"Points calculés: {len(ts)}")
    print(f"Altitude max: {max(hs):.3f} m")
    print(f"Vitesse max: {max(vs):.3f} m/s")
    print(f"Temps final: {ts[-1]:.3f} s")
    print(f"Volume eau final: {Vws[-1]:.6f} m³")
    
    # Afficher les premiers points
    print(f"\nPremiers points:")
    for i in range(min(5, len(ts))):
        print(f"t={ts[i]:.3f}s, h={hs[i]:.3f}m, v={vs[i]:.3f}m/s")
else:
    print("ERREUR: Aucun point calculé!")

# Test simple de txtGraph avec des données de test
print("\n=== TEST txtGraph ===")
test_t = [0, 1, 2, 3]
test_h = [0, 10, 20, 15]
txtGraph(test_t, test_h)
print("Fichier Output.txt créé")