import math

# ============================================================
# CONSTANTES PHYSIQUES ET PARAMÈTRES DU ROCKET À EAU
# ============================================================

g = 9.8               # Accélération de la gravité (m/s²)
rho_w = 1.0e3         # Densité de l'eau (kg/m³)
rho_atm = 1.23        # Densité de l'air (kg/m³)
Cd = 0.35             # Coefficient de traînée aérodynamique
gamma = 1.4           # Coefficient adiabatique de l'air
patm = 101325         # Pression atmosphérique (Pa)

# Paramètres géométriques et initiaux du rocket
#les noms avec o sont des paramètres à t = 0s (on peut pas mettre 0 dans la variable donc o)
# les paramètres qu on ne connait pas son notés float pour l'instant
D = float()               # Diamètre de la fusée (m)
De = float()            # Diamètre de la buse (m)
A = math.pi * (D / 2)**2     # Aire frontale du rocket (m²)
Ae = math.pi * (De / 2)**2   # Aire de la buse (m²)
V = float()          # Volume total du rocket (m³)
p_ino = float()           # Pression initiale à l'intérieur (Pa)
mb = float()             # Masse structurelle (kg)
Vwo = float()              # Volume d'eau initial dans la fusée
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
    Vao = V - Vwo
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
    F = -dmdt*v_e

    # --- Forces aérodynamiques
    FD, W = FD_W(v)

    # --- Accélération du rocket
    dvdt = (F - FD - W) / (mb + mw)

    # --- Variation de la pression interne
    dp_indt = -gamma * p_ino**((1+gamma)/gamma) * v_e * Ae / ((1-k) * V * p_ino**((1+gamma)/gamma))

    return F, dvdt, dp_indt
