import math
import matplotlib.pyplot as plt

# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.81               # Accélération de la gravité (m/s²)
rho_w = 1000.0         # Densité de l'eau (kg/m³)
rho_atm = 1.225        # Densité de l'air (kg/m³) à 15°C
Cd = 0.85              # Coefficient de traînée aérodynamique
gamma = 1.4            # Coefficient adiabatique de l'air
patm = 101325          # Pression atmosphérique (Pa)
R_air = 287.05         # Constante des gaz parfaits (J/kg·K)

# Paramètres géométriques
D = 0.1                # Diamètre de la fusée (m)
De = 0.008             # Diamètre de la buse (m)
A = math.pi * (D / 2)**2      # Aire frontale (m²)
Ae = math.pi * (De / 2)**2    # Aire de la buse (m²)
V = 0.0015             # Volume total (m³) = 1.5 L
p_ino = 500000         # Pression initiale (Pa) = 4 bars absolus
T_initial = 293.15     # Température initiale (20°C)
mb = 0.1               # Masse structurelle (kg)
Vwo = 0.0005           # Volume d'eau initial (m³) = 0.5 L

# Calculs dérivés
Vao = V - Vwo          # Volume d'air initial (m³)
m_air_initial = (p_ino * Vao) / (R_air * T_initial)  # Masse d'air initiale

# ============================================================
#  FONCTIONS CORRIGÉES
# ============================================================

def Weight(m):
    """Poids de la fusée (m = masse totale)"""
    return m * g

def Drag(v):
    """Force de traînée"""
    return 0.5 * rho_atm * abs(v) * v * Cd * A

def water_thrust(v_e):
    """Poussée pour l'EAU"""
    return rho_w * Ae * v_e**2

def air_thrust(m_dot, v_e):
    """Poussée pour l'AIR (formule générale)"""
    return m_dot * v_e

def internal_pressure_water(Vw):
    """Pression pendant la phase eau"""
    Va = V - Vw
    Vao = V - Vwo
    if Va > 0:
        return p_ino * (Vao / Va)**gamma
    return patm

def water_exit_velocity(p_in):
    """Vitesse d'éjection de l'EAU (Bernoulli simplifié)"""
    if p_in > patm:
        return math.sqrt(2 * (p_in - patm) / rho_w)
    return 0

def air_exit_velocity_realistic(p_in, T_air):
    """
    Vitesse d'éjection RÉALISTE de l'air
    Valeurs attendues: 100-250 m/s max
    """
    if p_in <= patm:
        return 0
    
    delta_p = p_in - patm
    PR = p_in / patm
    
    # Vitesse du son à la température de sortie
    T_exit = T_air * (patm / p_in) ** ((gamma - 1) / gamma)
    sound_speed = math.sqrt(gamma * R_air * T_exit)
    
    # Écoulement critique ?
    PR_critical = ((gamma + 1) / 2) ** (gamma / (gamma - 1))  # ≈ 1.893
    
    if PR >= PR_critical:
        # Écoulement critique - vitesse limitée au sonique
        v_e = sound_speed * 0.9  # ~90% du sonique
    else:
        # Écoulement sous-critique
        rho_exit = patm / (R_air * T_exit)
        v_e = math.sqrt(2 * delta_p / rho_exit) * 0.6  # Facteur de correction
    
    # Limite supérieure réaliste
    return min(v_e, 300)  # Jamais plus de 300 m/s

def air_mass_flow_rate_realistic(p_in, T_air):
    """
    Débit massique RÉALISTE d'air
    Valeurs attendues: 0.001 à 0.01 kg/s (1-10 g/s)
    """
    if p_in <= patm:
        return 0
    
    # Formule simplifiée mais réaliste
    delta_p = p_in - patm
    
    # Pour une buse de 8mm (Ae ≈ 5e-5 m²)
    # Coeff de débit typique pour l'air: Cd ≈ 0.6-0.9
    Cd_nozzle = 0.7
    
    # Densité à la sortie
    T_exit = T_air * (patm / p_in) ** ((gamma - 1) / gamma)
    rho_exit = patm / (R_air * T_exit)
    
    # Débit massique
    m_dot = Cd_nozzle * Ae * math.sqrt(2 * rho_exit * delta_p)
    
    return m_dot

# ============================================================
#  ÉQUATIONS CORRIGÉES
# ============================================================

def equation_vel_water(v, Vw, p_in):
    """Accélération pendant la phase EAU"""
    mw_current = rho_w * Vw
    m_total = mb + mw_current
    
    # Forces
    F_drag = Drag(v)
    F_weight = Weight(m_total)
    
    if Vw > 1e-6 and p_in > patm:
        v_e = water_exit_velocity(p_in)
        F_thrust = water_thrust(v_e)
    else:
        F_thrust = 0
    
    dv_dt = (F_thrust - F_drag - F_weight) / m_total
    return dv_dt

def equation_vel_air(v, p_in):
    """Accélération pendant la phase AIR"""
    # En phase air, plus d'eau
    m_total = mb
    
    # Forces résistives
    F_drag = Drag(v)
    F_weight = Weight(m_total)
    
    if p_in > patm:
        # Température de l'air
        T_air = T_initial * (p_in / p_ino) ** ((gamma - 1) / gamma)
        
        # Calculs réalistes
        v_e = air_exit_velocity_realistic(p_in, T_air)
        m_dot = air_mass_flow_rate_realistic(p_in, T_air)
        F_thrust = air_thrust(m_dot, v_e)
        
        # DEBUG - Afficher pour vérifier
        print(f"  Air: p={p_in/1e5:.2f}bar, T={T_air:.0f}K, "
              f"v_e={v_e:.0f}m/s, m_dot={m_dot*1000:.2f}g/s, "
              f"F={F_thrust:.2f}N")
    else:
        F_thrust = 0
    
    dv_dt = (F_thrust - F_drag - F_weight) / m_total
    return dv_dt

def equation_water_volume(v, Vw):
    """Variation du volume d'eau"""
    if Vw <= 1e-6:
        return 0
    
    p_in = internal_pressure_water(Vw)
    v_e = water_exit_velocity(p_in)
    
    if v_e > 0:
        return -Ae * v_e
    return 0

# ============================================================
#  RUNGE-KUTTA CORRIGÉ
# ============================================================

def Rungekutta(stepNbr=5000):
    """Simulation avec gestion correcte des phases"""
    
    # Variables pour suivre la phase
    phase = 'water'
    transition_time = 0
    
    def systeme_complet(t, y, current_phase):
        """
        y = [h, v, Vw, m_air, p_in]
        """
        h, v, Vw, m_air, p_in = y
        
        # Position
        dh_dt = v
        
        if current_phase == 'water':
            # PHASE EAU
            if Vw > 1e-6:
                # Calcul pression
                p_in = internal_pressure_water(Vw)
                dp_in_dt = 0  # La pression est calculée, pas intégrée
                
                # Volume d'eau
                dVw_dt = equation_water_volume(v, Vw)
                
                # Masse d'air constante
                dm_air_dt = 0
                
                # Vitesse
                dv_dt = equation_vel_water(v, Vw, p_in)
            else:
                # Transition vers phase air
                dVw_dt = 0
                dm_air_dt = 0
                dp_in_dt = 0
                dv_dt = 0
                
        else:
            # PHASE AIR
            dVw_dt = 0
            Vw = 0
            
            # Calcul de la pression (loi des gaz parfaits)
            Va = V  # Tout le volume disponible
            if m_air > 0:
                T_air = T_initial * (m_air / m_air_initial) ** ((gamma - 1) / gamma)
                p_in = (m_air * R_air * T_air) / Va
                p_in = max(p_in, patm)
                
                # Débit massique d'air
                m_dot = air_mass_flow_rate_realistic(p_in, T_air)
                dm_air_dt = -m_dot
                
                # Dérivée de pression (approximative)
                dp_in_dt = -m_dot * R_air * T_air / Va
            else:
                p_in = patm
                dm_air_dt = 0
                dp_in_dt = 0
            
            # Vitesse
            dv_dt = equation_vel_air(v, p_in)
        
        return [dh_dt, dv_dt, dVw_dt, dm_air_dt, dp_in_dt]
    
    # Conditions initiales
    y0 = [0.0, 0.0, Vwo, m_air_initial, p_ino]
    ts = [0.0]
    ys = [y0]
    
    deltaT = 0.001  # Pas plus petit pour précision
    i = 0
    
    print("=== DÉBUT SIMULATION ===")
    print(f"Pression initiale: {p_ino/1e5:.1f} bar")
    print(f"Volume eau: {Vwo*1000:.0f} mL")
    print(f"Masse air initiale: {m_air_initial*1000:.2f} g")
    print("-" * 50)
    
    while i < stepNbr and ys[-1][0] >= -0.1:
        i += 1
        t_current = ts[-1]
        y_current = ys[-1]
        
        # Détection de transition
        if phase == 'water' and y_current[2] <= 1e-6:
            phase = 'air'
            transition_time = t_current
            print(f"\n=== TRANSITION EAU → AIR à t={t_current:.3f}s ===")
            print(f"Pression au début phase air: {y_current[4]/1e5:.2f} bar")
        
        # RK4
        k1 = systeme_complet(t_current, y_current, phase)
        
        y_temp = [y_current[j] + deltaT/2 * k1[j] for j in range(5)]
        k2 = systeme_complet(t_current + deltaT/2, y_temp, phase)
        
        y_temp = [y_current[j] + deltaT/2 * k2[j] for j in range(5)]
        k3 = systeme_complet(t_current + deltaT/2, y_temp, phase)
        
        y_temp = [y_current[j] + deltaT * k3[j] for j in range(5)]
        k4 = systeme_complet(t_current + deltaT, y_temp, phase)
        
        # Mise à jour
        y_new = [y_current[j] + deltaT/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) 
                for j in range(5)]
        
        # Contraintes physiques
        y_new[2] = max(0, y_new[2])  # Vw >= 0
        y_new[3] = max(0, y_new[3])  # m_air >= 0
        y_new[4] = max(patm, y_new[4])  # p_in >= patm
        
        ts.append(t_current + deltaT)
        ys.append(y_new)
        
        # Critère d'arrêt pour phase air
        if phase == 'air' and y_new[4] <= patm * 1.01:
            print(f"Propulsion air terminée à t={t_current+deltaT:.3f}s")
            break
    
    # Résultats
    print(f"\n=== RÉSULTATS FINAUX ===")
    print(f"Durée totale: {ts[-1]:.3f} s")
    print(f"Hauteur max: {max([y[0] for y in ys]):.2f} m")
    print(f"Vitesse max: {max([y[1] for y in ys]):.2f} m/s")
    print(f"Pression finale: {ys[-1][4]/1e5:.3f} bar")
    
    # Graphiques
    plot_results(ts, ys)
    
    return ts, ys

def plot_results(times, states):
    """Affichage des résultats"""
    
    h = [s[0] for s in states]
    v = [s[1] for s in states]
    Vw = [s[2] * 1000 for s in states]  # mL
    m_air = [s[3] * 1000 for s in states]  # g
    p_in = [s[4] / 1000 for s in states]  # kPa
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    
    # Trajectoire
    axs[0, 0].plot(times, h, 'b-', linewidth=2)
    axs[0, 0].set_xlabel('Temps (s)')
    axs[0, 0].set_ylabel('Hauteur (m)')
    axs[0, 0].grid(True)
    axs[0, 0].set_title('Trajectoire')
    
    # Vitesse
    axs[0, 1].plot(times, v, 'r-', linewidth=2)
    axs[0, 1].set_xlabel('Temps (s)')
    axs[0, 1].set_ylabel('Vitesse (m/s)')
    axs[0, 1].grid(True)
    axs[0, 1].set_title('Vitesse')
    
    # Pression
    axs[1, 0].plot(times, p_in, 'orange', linewidth=2)
    axs[1, 0].axhline(y=patm/1000, color='r', linestyle='--', label='Pression atm.')
    axs[1, 0].set_xlabel('Temps (s)')
    axs[1, 0].set_ylabel('Pression (kPa)')
    axs[1, 0].grid(True)
    axs[1, 0].legend()
    axs[1, 0].set_title('Pression interne')
    
    # Volume eau & Masse air
    ax2 = axs[1, 0].twinx()
    ax2.plot(times, Vw, 'g-', linewidth=2, alpha=0.5)
    ax2.set_ylabel('Volume eau (mL)', color='g')
    ax2.tick_params(axis='y', labelcolor='g')
    
    # Masse d'air
    axs[1, 1].plot(times, m_air, 'purple', linewidth=2)
    axs[1, 1].set_xlabel('Temps (s)')
    axs[1, 1].set_ylabel('Masse d\'air (g)')
    axs[1, 1].grid(True)
    axs[1, 1].set_title('Masse d\'air restante')
    
    plt.suptitle('Simulation de fusée à eau - Phases eau et air')
    plt.tight_layout()
    plt.show()

# ============================================================
#  EXÉCUTION
# ============================================================

if __name__ == "__main__":
    print("Test des valeurs réalistes:")
    print("-" * 50)
    
    # Test de la poussée eau
    p_test = 400000  # 4 bars
    v_e_water = water_exit_velocity(p_test)
    F_water = water_thrust(v_e_water)
    print(f"EAU à {p_test/1e5:.1f} bar:")
    print(f"  v_e = {v_e_water:.1f} m/s")
    print(f"  F = {F_water:.1f} N")
    print()
    
    # Test de la poussée air
    T_test = 293
    v_e_air = air_exit_velocity_realistic(p_test, T_test)
    m_dot_air = air_mass_flow_rate_realistic(p_test, T_test)
    F_air = air_thrust(m_dot_air, v_e_air)
    print(f"AIR à {p_test/1e5:.1f} bar:")
    print(f"  v_e = {v_e_air:.1f} m/s")
    print(f"  m_dot = {m_dot_air*1000:.2f} g/s")
    print(f"  F = {F_air:.1f} N")
    print(f"  Rapport F_air/F_eau = {F_air/F_water*100:.1f}%")
    print("-" * 50)
    
    # Lancer la simulation
    ts, ys = Rungekutta(10000)