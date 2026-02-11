import math
import matplotlib.pyplot as plt

# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.81               # Accélération de la gravité (m/s²)
rho_w = 1.0e3         # Densité de l'eau (kg/m³)
rho_atm = 1.23        # Densité de l'air (kg/m³)
Cd = 0.5              # Coefficient de traînée aérodynamique (réduit pour fusée lisse)
patm = 101325         # Pression atmosphérique (Pa)

# Paramètres géométriques et initiaux
D = 0.1               # Diamètre de la fusée (m)
De = 0.008            # Diamètre de la buse (m)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
mb = 0.1              # Masse structurelle (kg)
Vwo = 0.0005          # Volume d'eau initial (m³) = 500 mL
rho_w = 1000          # Densité de l'eau (kg/m³)
mw0 = rho_w * Vwo     # Masse d'eau initiale

# ============================================================
# DONNÉES EXPÉRIMENTALES DE POUSSÉE 
# ============================================================

poussee_mesuree = [
    13.286, 15.31, 16.682, 16.975, 16.486, 13.482, 14.82, 16.616, 16.29,
    14.559, 14.07, 13.384, 14.984, 14.723, 12.568, 12.274, 14.363, 15.408,
    14.429, 15.408, 17.237, 15.604, 13.417, 13.417, 13.841, 10.968,
    8.323, 8.225, 8.16, 6.429, 6.625, 6.397, 6.886, 5.809, 6.331,
    5.613, 3.556, 4.111, 4.177, 2.936, 2.577, 2.773, 2.609,
    0, 0, 0, 0, 0, 0, 0  
]

# moving average to get rid of the noise
poussee_lissee = []
for i in range(len(poussee_mesuree)):
    if i == 0:
        poussee_lissee.append(poussee_mesuree[i])
    elif i == len(poussee_mesuree) - 1:
        poussee_lissee.append(poussee_mesuree[i])
    else:
        # Moyenne sur 3 points
        lisse = (poussee_mesuree[i-1] + poussee_mesuree[i] + poussee_mesuree[i+1]) / 3
        poussee_lissee.append(lisse)

# Use the new data
poussee_mesuree = poussee_lissee

# Calculate the total lenght of the experimental simulation
duree_propulsion_mesuree = len(poussee_mesuree) * 0.01  # secondes

print(f"Durée de propulsion: {duree_propulsion_mesuree:.2f} s")
print(f"Poussée moyenne pendant propulsion: {sum(poussee_mesuree[:43])/43:.2f} N")
print(f"Poussée maximale: {max(poussee_mesuree):.2f} N")


# Calculate the total impulsion caused by the thrust
impulsion_totale = sum([f * 0.01 for f in poussee_mesuree])
print(f"Impulsion totale: {impulsion_totale:.2f} N·s")

def get_thrust_from_measured_data(t):
    """
    Retourne la poussée mesurée au temps t.
    Les données sont disponibles pour t < duree_propulsion_mesuree.
    Pour t >= duree_propulsion_mesuree, retourne 0.
    """
    if t < 0:
        return 0
    
    # calculate the index of the measure because wwe took 1 measure every 0.01 s
    index = int(round(t / 0.01))
    
    if index < 0:
        return poussee_mesuree[0]
    elif index >= len(poussee_mesuree):
        return 0
    else:
        return poussee_mesuree[index]

# ============================================================
#  FORCES EN PRESENCE
# ============================================================

def Weight(mw):
    """
    mw : masse d'eau actuelle (kg)
    Calcule le poids de la fusée
    """
    return (mb + mw) * g

def Drag(v):
    """
    v : vitesse du rocket (m/s)
    Calcule la force de traînée
    """
    # La traînée est proportionnelle au carré de la vitesse
    if v > 0:
        return 0.5 * rho_atm * v**2 * Cd * A
    else:
        return -0.5 * rho_atm * v**2 * Cd * A  # Pour la descente

def equation_vel(v, mw, t):
    """
    Équation de vitesse utilisant les données mesurées de poussée
    """
    # Masse totale actuelle
    masse_totale = mb + mw
    
    # Forces en présence
    F_thrust = get_thrust_from_measured_data(t)
    F_drag = Drag(v)
    F_weight = Weight(mw)
    
    # Accélération (F = ma)
    if masse_totale > 0:
        dv_dt = (F_thrust - F_drag - F_weight) / masse_totale
    else:
        dv_dt = 0
    
    return dv_dt

def estimate_water_consumption_rate(t):
    """
    Estimate the water consumption rate based on the thrust data.
    Using this equation: F = ṁ * v_e où v_e ~ sqrt(2*ΔP/ρ)
    """
    F_thrust = get_thrust_from_measured_data(t)
    
    if F_thrust <= 0:
        return 0
    
    
    v_e_estimated = 30  # m/s - standard value
    
    
    m_dot = F_thrust / v_e_estimated
    
    # Volumic comsuption (kg/s -> m³/s)
    #V_dot = m_dot / rho_w
    
    return -m_dot 

# ============================================================
# MÉTHODE DE RUNGE-KUTTA 4ÈME ORDRE
# ============================================================

def RungeKutta_water_rocket(max_time=10.0, dt=0.01):
    """
    Simulation complète de la fusée à eau avec données mesurées
    """
    # Conditions initiales
    t = 0.0
    h = 0.0      # hauteur initiale
    v = 0.0      # vitesse initiale
    mw = mw0     # masse d'eau initiale
    
    # Listes pour stocker les résultats
    times = [t]
    heights = [h]
    velocities = [v]
    water_masses = [mw]
    thrusts = [get_thrust_from_measured_data(t)]
    accelerations = [0]


    def systeme_complet(t, y):
        h, v, mw = y
        dh_dt = v
        dv_dt = equation_vel(v, mw, t)
        dm_dt = estimate_water_consumption_rate(t)
        return [dh_dt, dv_dt, dm_dt]
    
    # Simulation
    while t < max_time and h >= -0.1:  # S'arrêter si on touche le sol
        # Stocker les valeurs actuelles
        current_thrust = get_thrust_from_measured_data(t)
        
        # État actuel
        y = [h, v, mw]
        

        




        k1 = systeme_complet(t, y)
        
        k2 = systeme_complet(t + dt/2, [y[j] + dt/2 * k1[j] for j in range(3)])
        
        k3 = systeme_complet(t + dt/2, [y[j] + dt/2 * k2[j] for j in range(3)])
        
        k4 = systeme_complet(t + dt, [y[j] + dt * k3[j] for j in range(3)])
        # Mise à jour
        y_new = [y[j] + dt/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) for j in range(3)]
            
        t += dt

        h_new = y_new[0]
        v_new = y_new[1]
        mw_new = max(0, mw-0.01)#y_new[2]


        
        # Stocker les nouvelles valeurs
        times.append(t)
        heights.append(h_new)
        velocities.append(v_new)
        water_masses.append(mw_new)
        thrusts.append(get_thrust_from_measured_data(t))
        
        # Mettre à jour pour l'itération suivante
        h, v, mw = h_new, v_new, mw_new
        
        # Arrêter si la fusée redescend et a déjà atteint une hauteur significative
        if t > 5 and v < -5 and h < max(heights)/10:
            break
    
    return times, heights, velocities, water_masses, thrusts

# ============================================================
# VISUALISATION
# ============================================================

def plot_results(times, heights, velocities, water_masses, thrusts):
    """
    Crée des graphiques complets des résultats
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Hauteur vs temps
    axes[0, 0].plot(times, heights, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Temps (s)')
    axes[0, 0].set_ylabel('Hauteur (m)')
    axes[0, 0].set_title('Trajectoire de la fusée')
    axes[0, 0].grid(True)
    axes[0, 0].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Marquer l'apogée
    max_h = max(heights)
    max_h_time = times[heights.index(max_h)]
    axes[0, 0].plot(max_h_time, max_h, 'ro', markersize=10)
    axes[0, 0].annotate(f'Apogée: {max_h:.1f} m', 
                       xy=(max_h_time, max_h),
                       xytext=(max_h_time + 0.5, max_h),
                       arrowprops=dict(arrowstyle='->'))
    
    # Vitesse vs temps
    axes[0, 1].plot(times, velocities, 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Temps (s)')
    axes[0, 1].set_ylabel('Vitesse (m/s)')
    axes[0, 1].set_title('Vitesse de la fusée')
    axes[0, 1].grid(True)
    axes[0, 1].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Marquer la vitesse maximale
    max_v = max(velocities)
    max_v_time = times[velocities.index(max_v)]
    axes[0, 1].plot(max_v_time, max_v, 'ro', markersize=10)
    axes[0, 1].annotate(f'Vmax: {max_v:.1f} m/s', 
                       xy=(max_v_time, max_v),
                       xytext=(max_v_time + 0.5, max_v),
                       arrowprops=dict(arrowstyle='->'))
    
    # Poussée vs temps
    axes[0, 2].plot(times[:len(thrusts)], thrusts, 'g-', linewidth=2)
    axes[0, 2].set_xlabel('Temps (s)')
    axes[0, 2].set_ylabel('Poussée (N)')
    axes[0, 2].set_title('Profil de poussée mesuré')
    axes[0, 2].grid(True)
    axes[0, 2].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Masse d'eau vs temps
    
    axes[1, 0].plot(times, water_masses, 'purple', linewidth=2)
    axes[1, 0].set_xlabel('Temps (s)')
    axes[1, 0].set_ylabel('Masse d\'eau (KG)')
    axes[1, 0].set_title('Consommation d\'eau')
    axes[1, 0].grid(True)
    axes[1, 0].set_ylim(bottom=0)
    
    # Accélération vs temps (dérivée de la vitesse)
    accelerations = []
    for i in range(1, len(velocities)):
        acc = (velocities[i] - velocities[i-1]) / (times[i] - times[i-1])
        accelerations.append(acc)
    
    # Ajuster la longueur
    acc_times = times[1:]
    if len(acc_times) > len(accelerations):
        acc_times = acc_times[:len(accelerations)]
    
    axes[1, 1].plot(acc_times, accelerations, 'orange', linewidth=2)
    axes[1, 1].set_xlabel('Temps (s)')
    axes[1, 1].set_ylabel('Accélération (m/s²)')
    axes[1, 1].set_title('Accélération de la fusée')
    axes[1, 1].grid(True)
    axes[1, 1].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[1, 1].axhline(y=-g, color='r', linestyle='--', alpha=0.5, label='-g')
    axes[1, 1].legend()
    
    # Diagramme phase (vitesse vs hauteur)
    axes[1, 2].plot(heights, velocities, 'b-', linewidth=2)
    axes[1, 2].set_xlabel('Hauteur (m)')
    axes[1, 2].set_ylabel('Vitesse (m/s)')
    axes[1, 2].set_title('Diagramme de phase')
    axes[1, 2].grid(True)
    axes[1, 2].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[1, 2].axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    # Marquer le début et la fin
    axes[1, 2].plot(heights[0], velocities[0], 'go', markersize=10, label='Départ')
    axes[1, 2].plot(heights[-1], velocities[-1], 'ro', markersize=10, label='Atterrissage')
    axes[1, 2].plot(max_h, 0, 'ko', markersize=8, label='Apogée')
    axes[1, 2].legend()
    
    plt.tight_layout()
    plt.show()
    
    return max_h, max_v, max_h_time

# ============================================================
# SIMULATION PRINCIPALE
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("SIMULATION DE FUSÉE À EAU AVEC POUSSÉE MESURÉE")
    print("=" * 60)
    
    # Paramètres de la simulation
    masse_totale_initiale = mb + mw0
    print(f"\nParamètres initiaux:")
    print(f"- Masse structure: {mb:.3f} kg")
    print(f"- Masse d'eau initiale: {mw0:.3f} kg")
    print(f"- Masse totale initiale: {masse_totale_initiale:.3f} kg")
    print(f"- Volume d'eau: {Vwo*1000:.0f} mL")
    print(f"- Poussée moyenne: {sum(poussee_mesuree)/len(poussee_mesuree):.2f} N")
    
    # Ratio poussée/poids initial
    F_moyenne = sum(poussee_mesuree[:30])/30  # Moyenne sur les 0.3 premières secondes
    ratio_poussee_poids = F_moyenne / (masse_totale_initiale * g)
    print(f"- Ratio poussée/poids initial: {ratio_poussee_poids:.2f}")
    
    # Lancer la simulation
    print("\nLancement de la simulation...")
    times, heights, velocities, water_masses, thrusts = RungeKutta_water_rocket(max_time=15.0, dt=0.01)
    
    # Analyser les résultats
    max_h, max_v, max_h_time = plot_results(times, heights, velocities, water_masses, thrusts)
    
    # Afficher les résultats
    print("\n" + "=" * 60)
    print("RÉSULTATS DE LA SIMULATION")
    print("=" * 60)
    print(f"Apogée atteinte: {max_h:.2f} m")
    print(f"Vitesse maximale: {max_v:.2f} m/s ({max_v*3.6:.1f} km/h)")
    print(f"Temps à l'apogée: {max_h_time:.2f} s")
    print(f"Durée totale de vol: {times[-1]:.2f} s")
    
    # Temps de propulsion effectif (poussée > 0)
    temps_propulsion = sum([1 for f in thrusts if f > 1]) * 0.01
    print(f"Temps de propulsion effectif: {temps_propulsion:.2f} s")
    
    # Calculer la vitesse moyenne pendant la propulsion
    v_propulsion = []
    for i, t in enumerate(times):
        if i < len(thrusts) and thrusts[i] > 1:
            v_propulsion.append(velocities[i])
    
    if v_propulsion:
        v_moy_propulsion = sum(v_propulsion) / len(v_propulsion)
        print(f"Vitesse moyenne pendant propulsion: {v_moy_propulsion:.2f} m/s")
    
    # Énergie et performances
    vitesse_atterrissage = abs(velocities[-1])
    print(f"Vitesse à l'atterrissage: {vitesse_atterrissage:.2f} m/s")
    
    # Temps caractéristiques
    print(f"\nTemps caractéristiques:")
    for i, (t, h, v) in enumerate(zip(times, heights, velocities)):
        if i % int(len(times)/10) == 0 and i > 0:
            print(f"  t={t:.1f}s: h={h:.1f}m, v={v:.1f}m/s")