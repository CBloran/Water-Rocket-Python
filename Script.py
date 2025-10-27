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

#========================================================
#######             FONCTION GRAPHIQUE           ########
#========================================================
def graphique_temporel(donnees_x, donnees_y, largeur=80, hauteur=20, couleur_points='rouge'):
    # Codes de couleurs ANSI
    couleurs = {
        'rouge': '\033[91m',
        'vert': '\033[92m',
        'jaune': '\033[93m',
        'bleu': '\033[94m',
        'magenta': '\033[95m',
        'cyan': '\033[96m',
        'blanc': '\033[97m',
        'reset': '\033[0m'
    }
    
    # Sélectionner la couleur
    couleur_code = couleurs.get(couleur_points, couleurs['rouge'])
    reset_code = couleurs['reset']
    
    if len(donnees_x) != len(donnees_y):
        print(f"Erreur: Les listes X ({len(donnees_x)}) et Y ({len(donnees_y)}) doivent avoir la même taille")
        return
    
    if len(donnees_x) == 0:
        print("Aucune donnée à afficher")
        return
    
    # Normaliser les données X et Y
    min_x = min(donnees_x)
    max_x = max(donnees_x)
    min_y = min(donnees_y)
    max_y = max(donnees_y)
    
    amplitude_x = max_x - min_x if max_x != min_x else 1
    amplitude_y = max_y - min_y if max_y != min_y else 1
    
    # Normaliser les positions dans la grille
    positions_x = [
        int((x - min_x) / amplitude_x * (largeur - 1))
        for x in donnees_x
    ]
    
    positions_y = [
        int((y - min_y) / amplitude_y * (hauteur - 3))  # -3 pour l'espace des légendes
        for y in donnees_y
    ]
    
    # Créer la grille
    grille = [[' ' for _ in range(largeur)] for _ in range(hauteur)]
    
    # Dessiner l'axe horizontal
    for i in range(largeur):
        grille[hauteur-2][i] = '─'
    
    # Dessiner l'axe vertical
    for i in range(hauteur-2):
        grille[i][0] = '│'
    
    # Origine
    grille[hauteur-2][0] = '└'
    
    # Dessiner la courbe (lignes en noir)
    for i in range(len(positions_x) - 1):
        x1 = positions_x[i]
        x2 = positions_x[i + 1]
        y1 = hauteur - 2 - positions_y[i]
        y2 = hauteur - 2 - positions_y[i + 1]
        
        # Dessiner une ligne entre les points consécutifs (sans couleur)
        if x1 != x2:  # Éviter la division par zéro
            for x in range(min(x1, x2), min(max(x1, x2) + 1, largeur)):
                t = (x - x1) / (x2 - x1)
                y = int(y1 + t * (y2 - y1))
                if 0 <= y < hauteur - 1 and x < largeur:
                    grille[y][x] = '•'
        
        # Marquer les points de données (AVEC COULEUR)
        if 0 <= y1 < hauteur - 1 and x1 < largeur:
            grille[y1][x1] = f"{couleur_code}●{reset_code}"
    
    # Dernier point (AVEC COULEUR)
    dernier_x = positions_x[-1]
    dernier_y = hauteur - 2 - positions_y[-1]
    if 0 <= dernier_y < hauteur - 1 and dernier_x < largeur:
        grille[dernier_y][dernier_x] = f"{couleur_code}●{reset_code}"
    
    # Afficher le graphique
    print("\n" + "GRAPHIQUE X-Y".center(largeur))
    print("─" * largeur)
    
    for i, ligne in enumerate(grille):
        # Convertir la liste en string pour l'affichage
        ligne_affichage = ''.join(ligne)
        if i == 0:
            print(f"{max_y:8.2f} │{ligne_affichage}")
        elif i == hauteur - 2:
            print(f"{min_y:8.2f} │{ligne_affichage}")
        else:
            print(f"        │{ligne_affichage}")
    
    print(" " * 9 + "└" + "─" * (largeur-1))
    print(f" " * 9 + f"{min_x:.2f}" + " " * (largeur-20) + f"{max_x:.2f}")
    
    # Statistiques
    print(f"\nPoints: {len(donnees_x)}")
    print(f"X: min={min_x:.2f}, max={max_x:.2f}")
    print(f"Y: min={min_y:.2f}, max={max_y:.2f}")
    print(f"Couleur des points: {couleur_points}")

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
    graphique_temporel(x_simple, y_simple, 40, 40, 'vert')
    return (ts, ys)


tsys = RungeKutta(y0, t0, stepsNbr, Tfinal)
