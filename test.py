import math

# Constantes physiques
patm = 101325  # Pression atmosphérique (Pa)
p_ino = 300000  # Pression initiale (Pa)
V = 0.001  # Volume total de la bouteille (m³)
Vwo = 0.0005  # Volume d'eau initial (m³)
rho_w = 1000  # Densité de l'eau (kg/m³)
mb = 0.1  # Masse de la fusée vide (kg)
Ae = 0.0001  # Aire de la section de sortie (m²)
A = 0.001  # Aire de la section transversale (m²)
gamma = 1.4  # Coefficient adiabatique pour l'air
g = 9.81  # Accélération gravitationnelle (m/s²)

# Calcul des volumes initiaux
Vao = V - Vwo  # Volume d'air initial
k0 = Vwo / V  # Ratio eau/volume initial

def FD_W(v, m_current=mb):
    """
    Calcule la traînée et le poids
    v: vitesse
    m_current: masse totale actuelle (kg)
    """
    # Coefficient de traînée simplifié
    Cd = 0.5
    rho_air = 1.225  # Densité de l'air (kg/m³)
    
    F_drag = 0.5 * Cd * rho_air * A * v**2
    # La traînée s'oppose toujours au mouvement
    if v != 0:
        F_drag = -abs(F_drag) * (v/abs(v)) if v != 0 else 0
    
    F_weight = m_current * g
    return F_drag, F_weight

def internal_pressure(Vw):
    """
    Vw : volume d'eau restant
    Calcul de la pression interne en fonction du volume d'eau restant.
    """
    if Vw < 0:
        Vw = 0
    
    Va = V - Vw  # volume d'air restant
    Vao = V - Vwo  # volume d'air initial

    if Va > 0 and Vao > 0:
        p_in = p_ino * (Vao / Va)**gamma  # loi adiabatique
        return max(p_in, patm)
    else:
        return patm

def water_exit_velocity(k, p_in):
    """
    k : ratio of remaining water volume to total volume
    p_in : internal pressure
    """
    try:
        if p_in <= patm:
            return 0
        
        # Formule corrigée
        denominator = rho_w * (1 - (Ae/A)**2)
        if denominator <= 0:
            return 0
            
        v_e = math.sqrt(2 * (p_in - patm) / rho_w)
        return v_e
    except:
        return 0

#=========================================================================
########           Méthode de Runge-Kutta complète           ##########
#=========================================================================

def RungeKutta_complet(t0, h0, v0, Vw0, dt=0.001, max_time=10):
    """
    Intégration complète avec pas de temps fixe
    """
    # Conditions initiales
    t = t0
    h = h0
    v = v0
    Vw = Vw0
    
    # Listes pour stocker les résultats
    t_list = [t]
    h_list = [h]
    v_list = [v]
    Vw_list = [Vw]
    phase_list = [1]  # 1: propulsion, 2: balistique
    
    # Boucle d'intégration
    phase = 1  # Commence en phase de propulsion
    i = 0
    
    while t < max_time and h >= 0 and i < 10000:
        i += 1
        
        # Déterminer la phase actuelle
        if Vw <= 0 and phase == 1:
            phase = 2  # Passage à la phase balistique
            print(f"Transition à la phase balistique à t={t:.3f}s, h={h:.2f}m, v={v:.2f}m/s")
        
        # Masse actuelle
        if phase == 1:
            m = mb + rho_w * max(Vw, 0)
        else:
            m = mb
        
        # Fonction du système d'équations
        def system(t, y):
            h_y, v_y, Vw_y = y
            
            # Calcul des forces
            if phase == 1 and Vw_y > 0:
                # Phase propulsion
                p_in = internal_pressure(Vw_y)
                k = Vw_y / V
                v_e = water_exit_velocity(k, p_in)
                F_thrust = rho_w * Ae * v_e**2 if v_e > 0 else 0
                F_drag, F_weight = FD_W(v_y, m)
                
                # Équations différentielles
                dh_dt = v_y
                dv_dt = (F_thrust - F_drag - F_weight) / m
                dVw_dt = -Ae * v_e if v_e > 0 else 0
                
            else:
                # Phase balistique
                F_drag, F_weight = FD_W(v_y, m)
                
                dh_dt = v_y
                dv_dt = (-F_drag - F_weight) / m
                dVw_dt = 0
            
            return [dh_dt, dv_dt, dVw_dt]
        
        # RK4
        y = [h, v, Vw]
        k1 = system(t, y)
        
        y_temp = [y[j] + 0.5*dt*k1[j] for j in range(3)]
        k2 = system(t + 0.5*dt, y_temp)
        
        y_temp = [y[j] + 0.5*dt*k2[j] for j in range(3)]
        k3 = system(t + 0.5*dt, y_temp)
        
        y_temp = [y[j] + dt*k3[j] for j in range(3)]
        k4 = system(t + dt, y_temp)
        
        # Mise à jour
        for j in range(3):
            y[j] += dt/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
        
        # Mettre à jour les variables
        t += dt
        h = max(y[0], 0)  # Éviter les hauteurs négatives
        v = y[1]
        Vw = max(y[2], 0)  # Éviter les volumes négatifs
        
        # Stocker les résultats
        t_list.append(t)
        h_list.append(h)
        v_list.append(v)
        Vw_list.append(Vw)
        phase_list.append(phase)
        
        # Conditions d'arrêt
        if h <= 0 and v < 0:
            print(f"Atterrissage à t={t:.3f}s")
            break
        
        if abs(v) < 0.01 and h > 0 and phase == 2:
            # Apogée approximative
            print(f"Apogée approximative à t={t:.3f}s, h={h:.2f}m")
    
    return t_list, h_list, v_list, Vw_list, phase_list

# Version alternative avec pas de volume adaptatif
def RungeKutta_volumeAdaptatif(t0, h0, v0, Vw0, max_steps=1000):
    """
    Version qui s'adapte au débit d'eau
    """
    # Conditions initiales
    t = t0
    h = h0
    v = v0
    Vw = Vw0
    
    t_list = [t]
    h_list = [h]
    v_list = [v]
    Vw_list = [Vw]
    
    phase = 1
    i = 0
    
    while i < max_steps and h >= 0:
        i += 1
        
        # Déterminer le pas de temps adaptatif
        if phase == 1 and Vw > 0:
            # En phase propulsion, adapter le pas au débit
            p_in = internal_pressure(Vw)
            k = Vw / V
            v_e = water_exit_velocity(k, p_in)
            
            if v_e > 0:
                # Pas basé sur le temps pour vider un petit volume
                dVw_target = Vwo / 1000  # Vider 0.1% du volume initial par pas
                dt = dVw_target / (Ae * v_e)
                dt = min(dt, 0.01)  # Limiter à 10ms max
            else:
                dt = 0.001
        else:
            # Phase balistique : pas fixe
            dt = 0.01
            phase = 2
        
        # Masse actuelle
        m = mb + rho_w * max(Vw, 0) if phase == 1 else mb
        
        # RK4
        def deriv(t, y):
            h_y, v_y, Vw_y = y
            
            if phase == 1 and Vw_y > 0:
                # Propulsion
                p_in = internal_pressure(Vw_y)
                k = Vw_y / V
                v_e = water_exit_velocity(k, p_in)
                F_thrust = rho_w * Ae * v_e**2 if v_e > 0 else 0
                F_drag, F_weight = FD_W(v_y, m)
                
                return [v_y, 
                       (F_thrust - F_drag - F_weight) / m,
                       -Ae * v_e if v_e > 0 else 0]
            else:
                # Balistique
                F_drag, F_weight = FD_W(v_y, m)
                return [v_y, (-F_drag - F_weight) / m, 0]
        
        y = [h, v, Vw]
        k1 = deriv(t, y)
        k2 = deriv(t + dt/2, [y[j] + dt/2*k1[j] for j in range(3)])
        k3 = deriv(t + dt/2, [y[j] + dt/2*k2[j] for j in range(3)])
        k4 = deriv(t + dt, [y[j] + dt*k3[j] for j in range(3)])
        
        # Mise à jour
        for j in range(3):
            y[j] += dt/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
        
        t += dt
        h = max(y[0], 0)
        v = y[1]
        Vw = max(y[2], 0)
        
        # Stockage
        t_list.append(t)
        h_list.append(h)
        v_list.append(v)
        Vw_list.append(Vw)
        
        # Transition de phase
        if phase == 1 and Vw <= 0:
            phase = 2
            print(f"Fin de l'eau à t={t:.3f}s")
        
        # Arrêt conditions
        if h <= 0 and v <= 0:
            break
    
    return t_list, h_list, v_list, Vw_list

#=========================================================================
########           Fonctions d'affichage           ##########
#=========================================================================

def txtGraph(xs: list, ys: list, FileName="Output.txt"):
    """
    Crée une représentation graphique textuelle
    """
    if not xs or not ys:
        print(f"Erreur: listes vides pour {FileName}")
        return
    
    with open(FileName, "w") as source_file:
        # Normalisation pour l'affichage
        if len(xs) != len(ys):
            min_len = min(len(xs), len(ys))
            xs = xs[:min_len]
            ys = ys[:min_len]
        
        # Créer une échelle verticale
        max_y = max(ys)
        min_y = min(ys)
        range_y = max_y - min_y
        
        if range_y > 0:
            scale = 50 / range_y  # 50 caractères de hauteur
        else:
            scale = 1
        
        for i in range(len(ys)):
            # Position verticale normalisée
            pos = int((ys[i] - min_y) * scale)
            line = " " * pos + "|" + f" t={xs[i]:.2f}s, y={ys[i]:.2f}m"
            source_file.write(line + "\n")
        
        # Informations résumées
        source_file.write("\n" + "="*60 + "\n")
        source_file.write(f"Durée: {xs[-1]:.2f} s | Hauteur max: {max(ys):.2f} m | Points: {len(xs)}\n")

def afficher_resultats(t_list, h_list, v_list, Vw_list):
    """
    Affiche un résumé des résultats
    """
    if not t_list:
        print("Aucun résultat à afficher")
        return
    
    print("\n" + "="*60)
    print("RÉSULTATS DE LA SIMULATION")
    print("="*60)
    print(f"Durée totale: {t_list[-1]:.3f} s")
    print(f"Hauteur maximale: {max(h_list):.3f} m")
    print(f"Vitesse maximale: {max(v_list):.3f} m/s")
    print(f"Vitesse à l'apogée: {v_list[h_list.index(max(h_list))]:.3f} m/s")
    print(f"Volume d'eau initial: {Vwo*1e6:.1f} ml")
    print(f"Temps de propulsion: {next((t for t, Vw in zip(t_list, Vw_list) if Vw <= 0), t_list[-1]):.3f} s")
    print(f"Nombre de points: {len(t_list)}")
    print("="*60)
    
    # Sauvegarde complète
    with open("simulation_complete.csv", "w") as f:
        f.write("t(s),h(m),v(m/s),Vw(m3)\n")
        for i in range(len(t_list)):
            f.write(f"{t_list[i]:.4f},{h_list[i]:.4f},{v_list[i]:.4f},{Vw_list[i]:.6f}\n")

#=========================================================================
########           Exécution principale           ##########
#=========================================================================

if __name__ == "__main__":
    print("Simulation de fusée à eau")
    print("="*50)
    
    # Méthode 1: Pas de temps fixe (plus stable)
    print("\nMéthode 1: Pas de temps fixe")
    t1, h1, v1, Vw1, phase1 = RungeKutta_complet(
        t0=0.0,
        h0=0.0,
        v0=0.0,
        Vw0=Vwo,
        dt=0.001,
        max_time=5
    )
    
    afficher_resultats(t1, h1, v1, Vw1)
    txtGraph(t1, h1, "hauteur_fixe.txt")
    txtGraph(t1, v1, "vitesse_fixe.txt")
    
    # Méthode 2: Pas adaptatif
    print("\nMéthode 2: Pas adaptatif")
    t2, h2, v2, Vw2 = RungeKutta_volumeAdaptatif(
        t0=0.0,
        h0=0.0,
        v0=0.0,
        Vw0=Vwo,
        max_steps=2000
    )
    
    afficher_resultats(t2, h2, v2, Vw2)
    txtGraph(t2, h2, "hauteur_adaptatif.txt")
    txtGraph(t2, v2, "vitesse_adaptatif.txt")
    
    print("\nSimulations terminées!")
    print("Fichiers générés:")
    print("- hauteur_fixe.txt, vitesse_fixe.txt")
    print("- hauteur_adaptatif.txt, vitesse_adaptatif.txt")
    print("- simulation_complete.csv")